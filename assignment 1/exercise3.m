close all;
clear all;
clc;

%% Exercise 3

% training data
floorArea = [50,52,24,150,80,76,77,120,115,65,61,30];
prices = [392000,245000,135000,1800000,579000,655000,653000,1276000,1108000,566000,477000,176000];
N = size(floorArea,2); % number of train samples

% % normaliza data
% floorArea = (floorArea - min(floorArea)) / (max(floorArea) - min(floorArea));
% prices = (prices - min(prices)) / (max(prices) - min(prices));

% plot data
figure;
scatter(floorArea, prices);
title('house price data');

%% linear model

F2 = [ones(size(floorArea,2),1), floorArea'];
Y = prices';

theta1 = inv(F2'*F2)*F2'*Y;
J1 = sum( (Y-F2*theta1).^2 ,1 ) / 2

% plot fit
figure;
floorsModel = (-0.1:0.001:1.1) * (max(floorArea)-min(floorArea)) + min(floorArea);
F1Model = [ones(size(floorsModel,2),1), floorsModel'];
pricesModel1 = F1Model*theta1;
plot(floorsModel,pricesModel1,'r');
hold on;
scatter(floorArea, prices);
title('linear model');

%% quadratic model

F2 = repmat( floorArea', [1 3] ) .^ (ones(N,1) * (0:2));
Y = prices';

theta2 = inv(F2'*F2)*F2'*Y;
J2 = sum( (Y-F2*theta2).^2 ,1 ) / 2

% plot fit
figure;
floorsModel = (-0.1:0.001:1.1) * (max(floorArea)-min(floorArea)) + min(floorArea);
F2Model = repmat( floorsModel', [1 3] ) .^ (ones(size(floorsModel')) * (0:2));
pricesModel2 = F2Model*theta2;
plot(floorsModel,pricesModel2,'r');
hold on;
scatter(floorArea, prices);
title('quadratic model');

%% polynomial model

degree = 5;
F3 = repmat( floorArea', [1 degree+1] ) .^ (ones(N,1) * (0:degree));
Y = prices';

theta3 = inv(F3'*F3)*F3'*Y;
J3 = sum( (Y-F3*theta3).^2 ,1 ) / 2

% plot fit
figure;
floorsModel = (-0.1:0.001:1.1) * (max(floorArea)-min(floorArea)) + min(floorArea);
F3Model = repmat( floorsModel', [1 degree+1] ) .^ (ones(size(floorsModel')) * (0:degree));
pricesModel3 = F3Model*theta3;
plot(floorsModel,pricesModel3,'r');
hold on;
scatter(floorArea, prices);
title('polynomial model');



