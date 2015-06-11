function [u] = superresolution_SIMONJENNI(g,D,lambda)
% input: g: double gray scaled image
%        D: downscaling matrix
% lambda: parameter % output: u: inpainted image

[MD, ND] = size(g);
[MND, MN] = size(D);
SRfactor = sqrt(MND/MN);
M = MD / SRfactor;
N = ND / SRfactor;
u = reshape(D'*g(:),M,N)/SRfactor^2;

% Fixing parameters for the primal-dual algorithm
tau = 1/(sqrt(lambda));
sigma = 1/(tau*8);
theta = 1;
num_iterations = 3000;

% Initialize variables
u_bar = u;
w = zeros([size(u),2]);
v = zeros(MD,ND);

tic;
for i=1:num_iterations
    
    % v-Step
    v = (v+sigma*(reshape(D*u_bar(:), MD, ND)-g))/(1+sigma/lambda);
    
    % w-Step
    tmp = w+sigma*grad(u_bar);
    norm = sqrt(sum(sum(sum(tmp(:,:,:).^2))));
    w = tmp/max(1,norm);

    % u-Step
    u_old = u;
    u = u+tau*div(w)-tau*reshape(D'*v(:),M,N);

    % u_bar update
    du = theta*(u-u_old);
    u_bar = u+theta*(u-u_old);
    
    if sum(du(:).^2) < 1e-7
        break;
    end
    
end
toc;
u = u_bar;

end

function grad_x = grad(x)
% Computes the gradient of x using forward differences and Neumann boundary
% conditions. 
grad_x = zeros([size(x) 2]);
grad_x(1:end-1,:,1) = x(2:end,:)-x(1:end-1,:);
grad_x(:,1:end-1,2) = x(:,2:end)-x(:,1:end-1);
end

function div_x = div(x)
% Computes the divergence of x using forward differences and Neumann 
% boundary conditions. It is chosen such that it is the adjoint to the
% gradient operator grad(x).
term1 = zeros(size(x,1), size(x,2));
term2 = term1;
term1(2:end-1,:,1) = x(2:end-1,:,1)-x(1:end-2,:,1);
term1(end,:,1) = x(end-1,:,1);
term2(:,2:end-1) = x(:,2:end-1,2)-x(:,1:end-2,2);
term2(:,end) = x(:,end-1,2);
div_x = term1+term2;
end

function cost = cost(x,y,z, D, g, lambda)
tmp = grad(x);
cost = (D*x(:) - g(:))'*y(:) - sum(sum(y.^2))/(2*lambda)+tmp(:)'*z(:);
end
