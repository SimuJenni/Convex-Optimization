close all;
clear all;
clc;

%% Exercise 1

% training data
floorArea = [50,52,24,150,80,76,77,120,115,65,61,30];
prices = [392000,245000,135000,1800000,579000,655000,653000,1276000,1108000,566000,477000,176000];
N = size(floorArea,2); % number of train samples

% normalize data
floorAreaN = (floorArea - min(floorArea)) / (max(floorArea) - min(floorArea));
pricesN = (prices - min(prices)) / (max(prices) - min(prices));

X = [ones(size(floorAreaN,2),1), floorAreaN'];
Y = pricesN';

% plot normalized data
figure;
scatter(floorArea, prices);
title('house price data');
 
%% batch gradient descent
%  
% t = 1;
% n = 100; % iterations
% theta = zeros(2,n);
% theta(:,1) = [0;0]; % init
% for i=2:n
%     err = Y-X*theta(:,i-1);
%     theta(:,i) = theta(:,i-1) + t * sum( repmat(err,[1,2]) .* X, 1)' / N;
% end
% thetaEnd = theta(:,n);
% 
% % plot error
% J = sum( (repmat(Y,[1,n])-X*theta).^2, 1 )/ 2;
% figure;
% plot(1:n,J);
% title('error vs iteration (batch gradient descent)');
% 
% % plot fit
% figure;
% floorsModel = 0:0.001:1;
% pricesModel = thetaEnd(1) + floorsModel*thetaEnd(2);
% plot(floorsModel,pricesModel,'r');
% hold on;
% scatter(floorAreaN, pricesN);
% title('house prices model (batch gradient descent)');
% 
% % plot theta succession
% figure;
% plot(theta(1,:),theta(2,:),'-*');

%% stochastic gradient descent

% t = 1;
% n = 200; % iterations
% theta = zeros(2,n);
% theta(:,1) = [0;0]; % init
% for i=2:n
%     idx = randi(N);
%     err = Y(idx)-X(idx,:)*theta(:,i-1);
%     theta(:,i) = theta(:,i-1) + t * err * X(idx,:)';
% end
% thetaEnd = theta(:,n);
% 
% % plot error
% J = sum( (repmat(Y,[1,n])-X*theta).^2, 1 )/ 2;
% figure;
% plot(1:n,J);
% title('error vs iteration (stochastic gradient descent)');
% 
% % plot fit
% figure;
% floorsModel = 0:0.001:1;
% pricesModel = thetaEnd(1) + floorsModel*thetaEnd(2);
% plot(floorsModel,pricesModel,'r');
% hold on;
% scatter(floorAreaN, pricesN);
% title('house prices model (stochastic gradient descent)');
% 
% % plot theta succession
% figure;
% plot(theta(1,end-100:end),theta(2,end-100:end),'-*');


%% Analytical solution

thetaA = inv(X'*X)*X'*Y;
J = sum( (Y-X*thetaA).^2 ,1 ) / 2;

% plot fit
figure;
floorsModel = 0:0.001:1;
pricesModel = thetaA(1) + floorsModel*thetaA(2);
plot(floorsModel,pricesModel,'r');
hold on;
scatter(floorAreaN, pricesN);
title('house prices model (analytical solution)');

% plot theta succession
figure;
plot(theta(1,end-100:end),theta(2,end-100:end),'-*');
hold on;
scatter(thetaA(1), thetaA(2),20,[1 0 0],'filled');
title('theta succession (stochastic gradient descent)');

