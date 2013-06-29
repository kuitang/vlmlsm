function [ gamma1,gamma2,gamma2_2,N1,N2,N2_2,seed,Am,Bm,theta,W,K,zeta, J, thisN1,thisN2,thisN2_2,L,U ] = gamma12( n,epsilon,Tmax,Wmax,edgeProb,assoc,seed,maxiter )
%[ gamma1,gamma2,gamma2_2,N1,N2,N2_2,seed,Am,Bm,theta,W,K,zeta, J, thisN1,thisN2,thisN2_2,L,U ] = gamma12( n,epsilon,Tmax,Wmax,edgeProb,assoc,seed,maxiter )
%   generates random model, computes gamma and numGraphPoints with method 1
%   using 1st deriv (new SODA method) and method 2 using 2nd deriv (old
%   AISTATS method). method 2_2 is improved 2nd deriv method for better a
%   not yet implemented repulsive edges for method 2

% allow non assoc for gamma2 - can do or not?
% return an array for various values of epsilon, maybe 1, 0.1, 0.01, 0.001
% allow random or fixed... fixed can use all W max, T max? or rand +/- max?

% assoc=1 means fully assoc model, 0 means general

FACeta=1; % FACeta >1 tries to speed up convergence by going further than 1* step

if nargin<8; maxiter=60; end
if nargin>=7
    rng(seed); 
else
    seed=-1;
end

fprintf('Running with n=%d, eps=%6.4f, Tmax=%d, Wmax=%d, p=%5.3f, assoc=%d, seed=%d\n',n,epsilon,Tmax,Wmax,edgeProb,assoc,seed);
%suggest T=0.4, Zmult=0.2, Jmult=
%zeta = randn(N,1) * 0.2; J = (rand(N,N)-0.5*GENERALNOTASSOC)*fac*3/8*10/N; J = J + J'; J = J - diag(diag(J));

% originally Mooij created what he called theta and J, then for nwbbp
% converted theta->zeta. I flipped so lower down, have my convention for
% theta, W
%seed=seed+1; rng(seed); fprintf('Running with seed of %d\n',seed);
zeta = (rand(n,1)-0.5) * 2*Tmax/2; % was randn which drew from std normal
J = rand(n,n);
if assoc~=1; J=J-0.5; J=2*J; end  % Wmax=2*Wmax?
%J = J*Wmax/8; J = J + J'; J = J - diag(diag(J)); % instead use J=triu(J,1) * Wmax/4; J=J+J'; ?
J=triu(J,1) * Wmax/4; J=J+J';

%J=ones(n,n)*Wmax/4; J=J-diag(diag(J)); % make all const Wmax ********
%zeta=zeros(n,1);

for i=1:n
    for j=i+1:n
        if rand>edgeProb % remove edge
            J(i,j)=0; J(j,i)=0;
        end
    end
end
[Am,Bm,allAm,allBm,timesm] = nmbbp(n,zeta,J,maxiter,FACeta);
[Aw,Bw,allAw,allBw,timesw,L,U] = nwbbp2(n,zeta,J,2,Am,Bm); % just run 1 iteration to get L, U
%bplot(timesm,allAm,allBm, timesw,allAw,allBw,seed,FACeta);

% So we have Ai, Bi bounds from Am, Bm
% We need to compute gamma2 which is based on 2nd deriv bound, i.e. AISTATS
% paper, and gamma1 which is based on 1st deriv new bound

theta = 2 * zeta - 2 * sum(J,2); % this gets to my theta from Mooij style zeta
W = 4 * J;
alpha = exp(abs(W)) - 1;

% Generate gamma=radius required, with 2 methods
% need points Ai + gamma, Ai + 3 gamma, ..., 1-Bi
% also generate numGraphPoints
%epsilon=0.1;

%%%%%%%%%%  Old method (AISTATS), constructs gamma2 using 2nd deriv bound
eta=min(Am,Bm);
K=zeros(n); 
z=sum(W~=0);
Sigma=sum(z)/n/n;
% For now assume all edges are assoc *** FIX!
% for i=1:n
%     for j=i+1:n
%         K(i,j)=eta(i)*eta(j)*(1-eta(i))*(1-eta(j))*(2*alpha(i,j)+1)/((alpha(i,j)+1)^2);
%     end
% end
a=0; b=0; a2=0; %a2 is new method using 2nd deriv from SODA
for i=1:n
    for j=i+1:n
        if W(i,j)~=0 % edge
            % old method
            if W(i,j)>0
                k=alpha(i,j)/(alpha(i,j)+1);
            else % repulsive
                k=-alpha(i,j);
            end
            K(i,j)=eta(i)*eta(j)*(1-eta(i))*(1-eta(j)) * (1-k^2);
            newa=k /K(i,j);
            a=max(a,newa);
            % new method
            Ai=Am(i); Bi=Bm(i);
            if W(i,j)>0
                newa2=alpha(i,j)/(1+alpha(i,j));
                Aj=Am(j); Bj=Bm(j);
            else % repulsive edge: flip j, now use result for assoc edge
                newa2=-alpha(i,j);
                Aj=Bm(j); Bj=Am(j);
            end
            k=newa2*newa2;
            if (1-Bi)<=Aj % i range <= j range
                fac=Bi*Aj-(1-Bi)*(1-Aj)*k;
            elseif (1-Bj)<=Ai % j range <= i range
                fac=Bj*Ai-(1-Bj)*(1-Ai)*k;
            elseif Ai<=Aj % overlap, i interval starts lower - BUT j RANGE COULD BE WITHIN i RANGE!
                poss1=Aj; poss2=min(1-Bi,1-Bj);
                fac=(1-k)*min( poss1*(1-poss1), poss2*(1-poss2) ); 
%                 mu=min([Aj,1-Bi, 1-Bj]); % middle elements; added the 3rd element
%                 oldfac=mu*(1-mu)*(1-k); if fac~=oldfac; disp('DIFFERENT'); end
%                 fprintf('newfac=%11.8f, oldfac=%11.8f\n',fac,oldfac);
            elseif Aj<=Ai % overlap, j interval starts lower
                poss1=Ai; poss2=min(1-Bi,1-Bj);
                fac=(1-k)*min( poss1*(1-poss1), poss2*(1-poss2) );
%                 mu=min([Ai,1-Bj, 1-Bi]); % middle elements
%                 oldfac=mu*(1-mu)*(1-k); if fac~=oldfac; disp('DIFFERENT'); end
%                 fprintf('newfac=%11.8f, oldfac=%11.8f\n',fac,oldfac);
            else
                disp('Error!')
                beep;
            end
            newa2=newa2/fac;
            a2=max(a2,newa2);
            
        end
    end
end
a=a/4;
for i=1:n
    newb=0;
    for j=1:n
        if W(i,j)~=0
            newb=newb+((alpha(i,j)+1)^2)/(2*alpha(i,j)+1);
        end
    end
    newb=(1-z(i)+newb) / ( eta(i) * (1-eta(i)) );
    b=max(b,newb);
end
a
a2
b
Omega=max(a,b); 
Lambda=n*Omega*sqrt(Sigma);
gamma2=sqrt(2*epsilon/Lambda/n);

%improved 2nd deriv method
Omega2=max(a2,b); Lambda2=n*Omega2*sqrt(Sigma); 
gamma2_2=sqrt(2*epsilon/Lambda2/n);

%%%%%%%%%%  New method (SODA), constructs gamma1 using 1st deriv bound
% This already works for general (not just fully assoc) models
D1bound=sqrt( sum( (sum(abs(W))'-log(L)-log(U)).^2 ) );
gamma1=epsilon/D1bound/sqrt(n);

N2=0; N1=0; N2_2=0;
thisN2=zeros(n,1); thisN1=zeros(n,1); thisN2_2=zeros(n,1);
for i=1:n
    frac=(1-Bm(i)-Am(i)-2*gamma2)/2/gamma2;
    thisN2(i)=round(frac+0.51)+1;
    frac=(1-Bm(i)-Am(i)-2*gamma2_2)/2/gamma2_2;
    thisN2_2(i)=round(frac+0.51)+1;
    frac=(1-Bm(i)-Am(i)-2*gamma1)/2/gamma1;
    thisN1(i)=round(frac+0.51)+1;
end
N2=sum(thisN2); N2_2=sum(thisN2_2); N1=sum(thisN1);

fprintf('Old second deriv method: gamma2=%8.6f, N2=%d\n',gamma2,N2);
fprintf('New first deriv method:  gamma1=%8.6f, N1=%d,    N2/N1=%d\n',gamma1,N1,N2/N1);
fprintf('Improved second deriv method: gamma2_2=%8.6f, N2_2=%d,    N2/N2_2=%d\n',gamma2_2,N2_2,N2/N2_2);
end
    
    
% Vary n, T, W, prob edge
        
        