clear all

load('data.mat');

% Resize blur 
k = rot90(f,2);
k = imresize(k,[5 5]);
k = k.*(k > 0);
k = k./sum(k(:));

% Resize sharp image
f = imresize(x,[51 51]);

% Synthetic blurred image 
g = conv2(f,k,'same');
[M, N] = size(f);
[MK, NK] = size(k);

%% Define matricies
lambda = 1e-4;
K = convmtx2(k,M,N);
Dx = convmtx2([1 -1],M,N);
Dy = convmtx2([1 -1]',M,N);

valid = zeros(M + MK - 1, N + NK - 1);
valid(floor(MK/2) + (1:M), floor(NK/2) + (1:M)) = 1;
K = K(logical(valid(:)),:);

valid = zeros(M + 1, N);
valid(1:M, :) = 1;
Dx = Dx(logical(valid(:)),:);

valid = zeros(M, N+1);
valid(:, 1:N) = 1;
Dy = Dy(logical(valid(:)),:);


%keyboard
%% Exact solution
fe = pinv(full(K'*K + lambda*(Dx'*Dx) + lambda*(Dy'*Dy)))*K'*g(:);
fe = reshape(fe,[M N]);
imagesc(fe)

%% Gradient descent
max_it = 5000;
it = 1;
alpha = 0.05;
beta = 0.05;
converged = false;
fg = g;
cost = @(x) sum((K*x(:) - g(:)).^2)/2 + lambda * sum(Dx*x(:).^2 + Dy*x(:).^2);
while it < max_it && ~converged
    Delta = K'*(K*fg(:)) + lambda*(Dx'*Dx + Dy'*Dy)*fg(:) - K'*g(:);
    
    % backtrack search
    t = 1;
    while(cost(fg(:) - t*Delta) > cost(fg(:)) - alpha*t*(Delta'*Delta))
        t = beta*t;
    end
    
    % Update
    fg = fg - t * reshape(Delta,M,N);
   
    % Stopping criterion
    if norm(Delta(:),2) < 3e-6*M*N
        converged = true;
    end
    
    it = it + 1;
    
    % visualization
    imagesc(fg)
    drawnow
end


%% Newton's method
max_it = 10;
it = 1;
converged = false;
alpha = 0.05;
beta = 0.05;
fn = g;
cost = @(x) sum((K*x(:) - g(:)).^2)/2 + lambda * sum(Dx*x(:).^2 + Dy*x(:).^2);
while it < max_it && ~converged
    Delta = K'*(K*fn(:)) + lambda*Dx'*(Dx*fn(:)) + lambda*Dy'*(Dy*fn(:)) - K'*g(:);
    Hess_inv = pinv(full(K'*K + lambda*(Dx'*Dx) + lambda*(Dy'*Dy)));
    
    % backtrack search
    t = 1;
    while(cost(fn(:) - t*Delta) > cost(fn(:)) - alpha*t*(Delta'*Delta))
        t = beta*t;
    end
    
    % Update
    fn = fn - t * reshape(Hess_inv * Delta,M,N);
    
    % Stopping criterion
    if (Delta' * Hess_inv * Delta)/2 < 3e-6*M*N
        converged = true;
    end
 
    it = it + 1;
    
    % visualization
    imagesc(fn)
    drawnow
end

