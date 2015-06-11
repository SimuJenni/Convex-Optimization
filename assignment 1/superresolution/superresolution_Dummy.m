function [u] = superresolution_Dummy(g,D,lambda)
% input: g: double gray scaled image
%        D: downscaling matrix
% lambda: parameter % output: u: inpainted image

% dummy output
[MD, ND] = size(g);
[MND, MN] = size(D);
SRfactor = sqrt(MND/MN);
M = MD / SRfactor;
N = ND / SRfactor;

% Initializing u including its boundaries
utmp = reshape(D'*g(:),M,N);
% u = zeros(M+2,N+2);
% u(2:end-1,2:end-1)=utmp;
u=padarray(utmp,[1 1], 'symmetric');

% Computing tau
bias = 1e-6;

% Gradient descent
max_it = 1000;
it = 1;
alpha = 0.1;
beta = 0.05;
t=0.001;
converged = false;
cost = @(x) lambda*sum(sum((D*x(:) - g(:)).^2))/2 + sum(sum(computeTau(x, bias)));
while it < max_it && ~converged
    tau = computeTau(u,bias);
    tau = reshape(tau,M+1,N+1);

    dTau1 = (2*u(2:end-1,2:end-1)-u(3:end,2:end-1)-u(2:end-1,3:end))./tau(2:end,2:end);
    dTau2 = (u(2:end-1,2:end-1)-u(1:end-2,2:end-1))./tau(1:end-1,2:end);
    dTau3 = (u(2:end-1,2:end-1)-u(2:end-1,1:end-2))./tau(2:end,1:end-1);
    deltaU = dTau1+dTau2+dTau3;
    
    Delta = lambda*D'*(D*reshape(u(2:end-1,2:end-1),[],1)-g(:))+deltaU(:);
    
    % backtrack search
    t = 1;
    while(cost(u(2:end-1,2:end-1) - t*reshape(Delta,M,N)) > cost(u(2:end-1,2:end-1)) - alpha*t*(Delta'*Delta))
        t = beta*t;
    end
    
    % Update
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) - t * reshape(Delta,M,N);
    u=padarray(u(2:end-1,2:end-1),[1 1], 'symmetric');
   
    % Stopping criterion
    if norm(Delta(:),2) < 3e-5*M*N
        converged = true;
    end
    
    it = it + 1;
    
    % visualization
    imagesc(u)
    drawnow
end
u=u(2:end-1,2:end-1);
end

function tau =computeTau(u, bias)
tau = sqrt((u(2:end,1:end-1)-u(1:end-1,1:end-1)).^2+(u(1:end-1,2:end)-u(1:end-1,1:end-1)).^2+bias);
end