function [A,B,allA,allB,times,L,U] = nwbbp(N,theta,J,maxiter,inpA,inpB)
% function [A,B,allA,allB,times] = nwbbp(N,theta,J,maxiter,inpA,inpB)
% new weller bbp, modified by AW 3/13
% Bethe bound propagation (http://arxiv.org/abs/1301.0015)
% 051013 added ability to pass in initial inpA and inpB (e.g. from MK)
% 052713 added output of L, U to work with 2nd deriv bound

% convert from (-1,+1) to (0,1) convention
zeta = 2 * theta - 2 * sum(J,2); % zeta is called theta in the paper
W = 4 * J;
alpha = exp(abs(W)) - 1;

allA=zeros(N,maxiter); allB=zeros(N,maxiter);
times=zeros(1,maxiter);
tstart=tic; % start timer

% initialize
Wsum = sum(W.*(W>0),2); V = sum(-W.*(W<0),2);
if nargin<5
    A = sigma(zeta-V);
else A=inpA;
end
if nargin<6
    B = 1 - sigma(zeta + Wsum);
else B=inpB;
end
allA(:,1)=A; allB(:,1)=B; 

times(1)=toc(tstart);

L=ones(N,1); U=ones(N,1); 
% main loop
for iter=1:maxiter-1 
  for i=1:N
    L_i = 1;
    U_i = 1;
    for j=1:N
      if W(i,j) > 0
        L_i = L_i * (1 + (alpha(i,j) * A(j)) / (1 + alpha(i,j) * (1-B(i)) * (1-A(j))));
        U_i = U_i * (1 + (alpha(i,j) * B(j)) / (1 + alpha(i,j) * (1-A(i)) * (1-B(j))));
      elseif W(i,j) < 0
        L_i = L_i * (1 + (alpha(i,j) * B(j)) / (1 + alpha(i,j) * (1-B(i)) * (1-B(j))));
        U_i = U_i * (1 + (alpha(i,j) * A(j)) / (1 + alpha(i,j) * (1-A(i)) * (1-A(j))));  
      end
    end
    A(i) = 1 / (1 + exp(-zeta(i) + V(i)) / L_i);
    B(i) = 1 / (1 + exp(zeta(i) + Wsum(i)) / U_i);
    L(i)=L_i; U(i)=U_i;
  end
  %fprintf('Finished an iteration and computed L etc.\n');
  allA(:,iter+1)=A; allB(:,iter+1)=B;
  times(iter+1)=toc(tstart);
end

%[A,1-B]

return

function y=sigma(x)
  y = 1 ./ (1 + exp(-x));
return
