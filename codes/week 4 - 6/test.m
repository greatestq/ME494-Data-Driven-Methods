clear all; close all; clc;

%% Generate Data for f(t) = t^2
tspan = linspace(0, 5, 500); % Time span
f = tspan.^2; % Define the function f(t) = t^2

%% Compute Derivative
dfdt = 2 * tspan; % Compute the derivative of f(t)

%% Build library and compute sparse regression
n = 1; % Number of states
polyorder = 3; % Maximum polynomial order for the library
Theta = poolData(f, n, polyorder); % Build the library
lambda = 0.025; % Sparsification parameter
Xi = sparsifyDynamics(Theta, dfdt, lambda, n); % Compute sparse regression

%% f(t) reconstruction
fdot_sindy = Theta*Xi;
f_sindy = cumtrapz(tspan, fdot_sindy);
%f_sindy2 = cumtrapz(fdot_sindy);
%% Plotting
figure;
subplot(1,2,1)
plot(tspan, dfdt, 'LineWidth', 2); hold on;
plot(tspan, Theta * Xi, '--', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Derivative of f(t)');
legend('True Derivative', 'SINDY Approximation');
title('SINDY Approximation for f(t) = t^2');
grid on;

subplot(1,2,2)
plot(tspan, f, 'LineWidth',2); hold on;
plot(tspan,f_sindy, '--', 'LineWidth',2);
xlabel('Time (t)');
ylabel('f(t)');
legend('True Function', 'SINDY Approximation');
grid on;
