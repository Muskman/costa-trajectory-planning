clc
close all
clear all

% generate environment with currents
[X,Y] = meshgrid(-2:0.2:2);
Z = X .* exp(-X.^2 - Y.^2);

[U,V] = gradient(Z,0.2,0.2);
figure; quiver(X,Y,U,V)

%%
T = 50;
dt = 0.1;
v = 0;                              % std of noise in disturbance perdictor

x_init = zeros(2*T,1);              % initial point

x = x_init; 
e = v*randn(2*T,1);                 % generate realization of noise 
[g, grad_f] = disturbance(x,e,dt);  % disturbance prediction and gradient of objective function


