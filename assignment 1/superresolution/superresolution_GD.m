function [u] = superresolution_SIMONJENNI(g,D,lambda)
% input: g: double gray scaled image
%        D: downscaling matrix
% lambda: parameter % output: u: inpainted image

[MD, ND] = size(g);
[MND, MN] = size(D);
SRfactor = sqrt(MND/MN);
M = MD / SRfactor;
N = ND / SRfactor;

% Initializing u
u = reshape(D'*g(:),M,N)/SRFactor^2;

% Bias term for computing tau
bias = 1e-6;

% Gradient descent parameters
max_it = 5000;
it = 1;
converged = false;

% Backtrack search parameters
alpha = 0.01;
beta = 0.05;
cost = @(x) lambda*sum(sum((D*x(:) - g(:)).^2))/2 + sum(sum(computeTau(x, bias)));
while it < max_it && ~converged
    % Computing tau
    tau = computeTau(u,bias);
    tau = reshape(tau,M+1,N+1);
    
    % Computation of gradient terms using mirror boundary assumptions on u
    dTau1 = (2*u-u([2:end end],:)-u(:,[2:end end]))./tau(2:end,2:end);
    dTau2 = (u-u([1 1:end-1],:))./tau(1:end-1,2:end);
    dTau3 = (u-u(:,[1 1:end-1]))./tau(2:end,1:end-1);
    deltaU = dTau1+dTau2+dTau3;
    
    Delta = lambda*D'*(D*reshape(u,[],1)-g(:))+deltaU(:);
    
    % backtrack search
    t = 1;
    while(cost(u - t*reshape(Delta,M,N)) > cost(u) - alpha*t*(Delta'*Delta))
        t = beta*t;
    end
    
    % Update u
    u = u - t * reshape(Delta,M,N);
   
    % Stopping criterion
    if norm(Delta(:),2) < 1e-6*M*N
        converged = true;
    end
    
    it = it + 1;
    
%     % visualization
%     imagesc(u)
%     drawnow
end
end

% Computes tau using mirrored boundaries for u
function tau =computeTau(u, bias)
tau = sqrt((u([1:end end],[1 1:end])-u([1 1:end],[1 1:end])).^2+...
    (u([1 1:end],[1:end end])-u([1 1:end],[1 1:end])).^2+bias);
end