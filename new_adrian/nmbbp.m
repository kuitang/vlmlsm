function [A,B,allA,allB,times] = nmbbp(N,theta,J,maxiter,stepscale)
% function [A,B,allA,allB,times] = nmbbp(N,theta,J,maxiter)
% new mooij bbp, modified by AW 3/13
% Improved bound calculation algorithm
% section IIIE of http://cs.ru.nl/~jorism/articles/04385778.pdf
% take step size * stepscale, see if improves

if nargin<5; stepscale=1; end % default to multiple of 1

all_eta_lo=zeros(N,N,maxiter); all_eta_hi=zeros(N,N,maxiter);
allA=zeros(N,maxiter); allB=zeros(N,maxiter);
times=zeros(1,maxiter);
tstart=tic; % start timer

% Initialize bounds on cavity fields
eta_lo = zeros(N,N);
eta_hi = zeros(N,N);
for i=1:N
  for j=1:N
    if J(i,j)
      eta_lo(i,j) = -Inf;
      eta_hi(i,j) = Inf;
    end
  end
end

% Propagate bounds on cavity fields
for iter=1:maxiter
  for i=1:N
    for j=1:N
      if J(i,j)
        eta_lo_new(i,j) = theta(i);
        eta_hi_new(i,j) = theta(i);
        for k=1:N
          if J(i,k) && k ~= j
            a = atanh(tanh(J(k,i)) * tanh(eta_lo(k,i)));
            b = atanh(tanh(J(k,i)) * tanh(eta_hi(k,i)));
            eta_lo_new(i,j) = eta_lo_new(i,j) + min([a b]);
            eta_hi_new(i,j) = eta_hi_new(i,j) + max([a b]);
          end
        end
      end
    end
  end
  if iter==1
      eta_lo = eta_lo_new;
      eta_hi = eta_hi_new;
  else
      change=stepscale*(eta_lo_new-eta_lo);
      eta_lo = eta_lo+change;
      change=stepscale*(eta_hi_new-eta_hi);
      eta_hi = eta_hi+change;
  end
  all_eta_lo(:,:,iter)=eta_lo; all_eta_hi(:,:,iter)=eta_hi;
  times(iter)=toc(tstart);
end

% convert each iteration. note we stopped the timer above already so giving
% this version the benefit of that
tic
for iter=1:maxiter
    
% Calculate beliefs in tanh parameterization
beta_lo = theta;
beta_hi = theta;
for i=1:N
  for k=1:N
    if J(i,k)
      a = atanh(tanh(J(k,i)) * tanh(all_eta_lo(k,i,iter)));
      b = atanh(tanh(J(k,i)) * tanh(all_eta_hi(k,i,iter)));
      beta_lo(i) = beta_lo(i) + min([a b]);
      beta_hi(i) = beta_hi(i) + max([a b]);
    end
  end
end

% Calculate A,B bounds on pseudomarginal
allA(:,iter) = (tanh(beta_lo) + 1) / 2;
allB(:,iter) = 1 - ((tanh(beta_hi) + 1) / 2);

end
A=allA(:,maxiter); B=allB(:,maxiter);
%[A,1-B]
fprintf('To handle the conversions at the end takes: '); toc

return
