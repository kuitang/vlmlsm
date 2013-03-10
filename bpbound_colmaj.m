function [A,B] = bpbound(N,theta,W,maxiter)
% function [A,B] = bpbound(N,theta,J,maxiter)
% Improved bound calculation algorithm
% section IIIE of http://cs.ru.nl/~jorism/articles/04385778.pdf

% Initialize bounds on cavity fields
eta_lo = zeros(N,N);
eta_hi = zeros(N,N);

% Convert from {0, 1} to {-1, 1}
% A bit more cumbersome than just inverting Mooij's equations; just write
% out the grid notation
J = 0.25 * W;
eta = 0.5*theta + sum(J, 2);

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
  for j=1:N
    for i=1:N
      if J(i,j)
        eta_lo_new(i,j) = eta(i);
        eta_hi_new(i,j) = eta(i);
        for k=1:N
          if J(k,i) && k ~= j
            a = atanh(tanh(J(k,i)) * tanh(eta_lo(k,i)));
            b = atanh(tanh(J(k,i)) * tanh(eta_hi(k,i)));
            eta_lo_new(i,j) = eta_lo_new(i,j) + min([a b]);
            eta_hi_new(i,j) = eta_hi_new(i,j) + max([a b]);
            
            if iter <= 3
                fprintf(1, 'i = %d, j = %d, k = %d\n', i - 1, j - 1, k - 1);
                fprintf(1, 'a = %g, b = %g, etaLo = %g, etaHi = %g, etaLoNew = %g, etaHiNew = %g\n', ...
                            full(a), full(b), eta_lo(k,i), eta_hi(k,i), eta_lo_new(i,j), eta_hi_new(i,j));                    
            end
            
          end
        end
      end
    end
  end
    
  eta_lo = eta_lo_new;
  eta_hi = eta_hi_new;
  
%  if iter < 3
%      [eta_lo(:) eta_hi(:) eta_lo_new(:) eta_hi_new(:)]
%  end
end
disp('eta')
[eta_lo(:) eta_hi(:)]

% Calculate beliefs in tanh parameterization
beta_lo = eta;
beta_hi = eta;
for i=1:N
  for k=1:N
    if J(i,k)
      fprintf(1, 'i = %d, j = %d\n', i, j);

      a = atanh(tanh(J(k,i)) * tanh(eta_lo(k,i)));
      b = atanh(tanh(J(k,i)) * tanh(eta_hi(k,i)));
      beta_lo(i) = beta_lo(i) + min([a b]);
      beta_hi(i) = beta_hi(i) + max([a b]);
    end
  end
end

disp('beta')
[beta_lo beta_hi]

% Calculate A,B bounds on pseudomarginal
A = (tanh(beta_lo) + 1) / 2;
B = 1 - ((tanh(beta_hi) + 1) / 2);

[A,1-B]

return
