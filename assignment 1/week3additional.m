clear all

%% Smoothing

% 1D example
N = 100; % number of samples
x = linspace(0,2,N); % function domain
f = sin(x*2*pi)+1e-1*randn(1,N); % function

figure(1)
clf
plot(f,'linewidth',4)
hold all
NL = 3;
lambdas = logspace(-2,0,NL);
legtext = cell(NL,1);
legtext{1} = 'f';

%fwd_diff = @(x) 

Nit = 1000; % number of iterations
for lam=1:NL
    lambda = lambdas(lam);
    u = f;
    e = 2e-1;
    for it = 1:Nit
        
        % periodic
        u_back = u([2:end 1]);
        u_forw = u([end 1:end-1]);
        u_back_2 = u([3:end 1 2]);
        u_forw_2 = u([end-1  end 1:end-2]);
        
        % Mirrored 
        % Equivalent to Neumann where u'[y] = 0
        % for y at the signal boundary 
%         u_back = u([2:end end]);
%         u_forw = u([1 1:end-1]);
%         u_back_2 = u([3:end end end-1]);
%         u_forw_2 = u([2  1 1:end-2]);

        % Dirichlet
%         ue = [0 0 u 0 0];
%         u_back = ue(2:end-3);
%         u_forw = ue(4:end-1);
%         u_back_2 = ue(1:end-4);
%         u_forw_2 = ue(5:end);

        % forward
        uxx = 2 * (2 * u - u_forw - u_back);
        % backward
        %uxx = 2 * (2 * u - u_forw - u_back);
        % central
        %uxx = 2 * u - u_forw_2 - u_back_2;
        
        u = u-e*(lambda*(u-f)+uxx);
        
    end
    plot(u,'linewidth',3)
    legtext{lam+1} = ['u with \lambda = ' num2str(lambda)];
end
hold off
xlabel('x','fontsize',20);
ylabel('f,u','fontsize',20);
legend(legtext);
set(gca,'fontsize',20);
drawnow

%% 1D TV

% 1D example
N = 100; % number of samples
x = linspace(0,2,N); % function domain
f = min(0.3,max(-0.3,sin(x*2*pi)))+1e-1*randn(1,N); % function

figure(1)
clf
plot(f,'linewidth',4)
hold all
NL = 3;
lambdas = logspace(-1,1,NL);
legtext = cell(NL,1);
legtext{1} = 'f';



Nit = 5000; % number of iterations
for lam=1:NL
    lambda = lambdas(lam);
    u = f;
    e = 1e-3;
    for it = 1:Nit
        
        
        
        % periodic
        u_back = u([2:end 1]);
        u_forw = u([end 1:end-1]);
        u_back_2 = u([3:end 1 2]);
        u_forw_2 = u([end-1  end 1:end-2]);
        
        % Mirrored 
        % Equivalent to Neumann where u'[y] = 0
        % for y at the signal boundary 
%         u_back = u([2:end end]);
%         u_forw = u([1 1:end-1]);
%         u_back_2 = u([3:end end end-1]);
%         u_forw_2 = u([2  1 1:end-2]);

        % Dirichlet
%         ue = [0 0 u 0 0];
%         u_back = ue(2:end-3);
%         u_forw = ue(4:end-1);
%         u_back_2 = ue(1:end-4);
%         u_forw_2 = ue(5:end);
        
        % forward
        uxx = -sign( u_back - u ) + sign(u - u_forw);
        % backward
        %uxx = sign( u - u_forw ) - sign(u_back - u);
        % central
        %uxx = sign( (u - u_forw_2)/2 ) - sign((u_back_2 - u)/2);
        
        u = u-e*(lambda*(u-f)+uxx);
    end
    
    plot(u,'linewidth',3)
    legtext{lam+1} = ['u with \lambda = ' num2str(lambda)];
end
hold off
xlabel('x','fontsize',20);
ylabel('f,u','fontsize',20);
legend(legtext);
set(gca,'fontsize',20);
drawnow