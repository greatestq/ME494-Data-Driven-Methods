% Clear workspace
clear; clc;

% directory where data is stored
DIR = './data';

paramsFile = fullfile(DIR,'airfoilDNS_parameters.h5');

dt_field = h5read(paramsFile,'/dt_field'); % timestep for field variables (velocity and vorticity)
dt_force = h5read(paramsFile,'/dt_force'); % timestep for scalar quantities
Re = h5read(paramsFile,'/Re');
FreqsAll = h5read(paramsFile,'/frequencies'); % pitching frequencies
alpha_p = h5read(paramsFile,'/alpha_p'); % pitching amplitude (deg)
alpha_0s = h5read(paramsFile,'/alpha_0s'); % base angles of attack (deg) (25 and 30)
pich_axis = h5read(paramsFile,'/pitch_axis'); % 0.5, midchord pitching

%h5disp(paramsFile)

%%
% load and plot snapshots for various airfoil kinematics
% each simulation contains 401 snapshots, with uniform timestep of 0.1 c/U_\infty time units
% load and plot snapshots
tstep = 200; % timestep to plot

% load spatial grid, and time vectors
filenameGrid = fullfile(DIR,'airfoilDNS_grid.h5');

x = h5read(filenameGrid,'/x');
y = h5read(filenameGrid,'/y');
nx = length(x);
ny = length(y);

figure

filename = fullfile(DIR,"airfoilDNS_a25f0p05.h5");
uy = h5read(filename,'/uy');

%Time varying airfoil coordinates
xa = h5read(filename,'/xa');
ya = h5read(filename,'/ya');
t_field = h5read(filename,'/t_field');
t_force = h5read(filename,'/t_force');
nt = length(t_field);

[~,hc] = contourf(x,y,squeeze(uy(:,:,tstep)).', linspace(-1.5,1.5,80));
set(hc,'LineStyle','none');
caxis([-0.5,0.5]);
hold on
plot(xa(:,tstep),ya(:,tstep),'k-')

colorbar;
axis image;
axis equal
title(['$u_y, \ \alpha_0 = 25 ^\circ$, $\alpha_P =  ',num2str(alpha_p) '^\circ$, $f_P = 0.05','$' ], 'interpreter', 'latex','fontsize',16)
xlabel('$x/c$','Interpreter','Latex');
ylabel('$y/c$','Interpreter','Latex');


%% Run POD Analysis

meanSub = 1;
filename = fullfile(DIR,"airfoilDNS_a25f0p05.h5");
ux = h5read(filename,'/ux');
uy = h5read(filename,'/uy');
uxreshape = reshape(ux,nx*ny,nt);
uyreshape = reshape(uy,nx*ny,nt);
data = [uxreshape; uyreshape];
if meanSub
    dataMean = mean(data,2);
    data = data-dataMean*ones(1,nt);
end

[U,S,V] = svd(data,'econ');

%% Sindy Practice
temporalAmplitude = V(:,6) * S(6,6);
f = temporalAmplitude.';

%Calculating the derivative of the temporal amplitude
%xDot = diff(x, 1, 1);
dt= mean(diff(t_field));
fdot = gradient(f, dt);

n = 6;
polyorder = 2;
Theta = poolData(f, n, polyorder);
lambda = 0.025;
Xi = sparsifyDynamics(Theta, fdot, lambda, n);

%% f(t) reconstruction
fdot_sindy = Theta*Xi;
f_sindy = cumtrapz(t_field.',fdot_sindy);
%f_sindy2 = cumtrapz(fdot_sindy);
figure;
subplot(1,3,1)
plot(t_field.', fdot, 'LineWidth', 2); hold on;
plot(t_field.', Theta * Xi, '--', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Derivative of f(t)');
legend('True Derivative', 'SINDY Approximation');
title('SINDY Approximation for first 6 temporal amplitudes');
grid on;

subplot(1,3,2)
plot(t_field.', f, 'LineWidth',2); hold on;
%plot(t_field.'-t_field(1), f_sindy, '--', 'LineWidth',2);
xlabel('Time (t)');
ylabel('f(t)');
legend('True Function');
grid on;

subplot(1,3,3)
plot(t_field.', f_sindy + f(1), '--', 'LineWidth',2, color='#D95319');
legend('SINDY Approximation')
grid on;

%% Individual temporal amplitudes
polyorder = 3;
lambda = 0.0135;
figure;
for i = 1:6
    n = 2;
    temporalAmplitude = V(:,i) * S(i,i);
    f = temporalAmplitude.';
    dt = mean(diff(t_field));
    fdot = gradient(f,dt);
    Theta = poolData(f, n, polyorder);
    Xi = sparsifyDynamics(Theta, fdot, lambda, n);

    fdot_sindy = Theta*Xi;
    f_sindy = cumtrapz(t_field.', fdot_sindy);

    subplot(2,6,i)
    plot(t_field.', fdot, 'LineWidth', 2); hold on;
    plot(t_field.', Theta * Xi, '--', 'LineWidth', 2);
    xlabel('Time (t)');
    ylabel('Derivative of f(t)');
    legend('True Derivative', 'SINDY Approximation');
    %title('SINDY Approximation for first', num2str(i), 'temporal amplitudes');
    grid on;
    title([num2str(i), 'th temporal amplitude']);  % Add the title
    xlim([50,90]);

    subplot(2,6,i+6)
    plot(t_field.', f, 'LineWidth',2); hold on;
    plot(t_field.', f_sindy + f(1), '--', 'LineWidth', 2, color='#D95319');
    xlabel('Time (t)');
    ylabel('f(t)');
    legend('True Function', 'SINDY derived');
    grid on;
    xlim([50,90]);
    ylim([-15,20])
end

%%
n = 6;
X = zeros(index, size(t_field,1));
Xdot = X;
dt = mean(diff(t_field));
for i = 1:n
    temporalAmplitude = V(:,i) * S(i,i);
    f = temporalAmplitude.';
    fdot = gradient(f,dt);
    X(i,:) = f;
    Xdot(i,:) = fdot;
end

polyorder = 2;
lambda = 0.0135;

Theta = poolData(X, n, polyorder);
Xi = sparsifyDynamics(Theta, Xdot, lambda, n);

fdot_sindy = Theta*Xi;

figure;
for j = 1:n
    f_sindy= cumtrapz(t_field.', fdot_sindy(j,:));
    subplot(2,6,j);
    plot(t_field.', fdot, 'LineWidth', 2); hold on
    plot(t_field.', fdot_sindy(j,:), '--k', 'LineWidth',2);
    xlabel('Time(t)');
    ylabel('Derivative of f(t)');
    legend('True Derivative', 'SINDY Approximation');
    grid on
    title([num2str(j), 'th temporal amplitude']);

    subplot(2,6,j+6)
    plot(t_field.', X(j,:), 'LineWidth', 2); hold on;
    plot(t_field.', f_sindy + X(j,1), '--k', 'LineWidth', 2);
    xlabel('Time(t)');
    ylabel('f(t)')
    legend('True Function', 'SINDY derived');
    grid on
end

%%
n = 6;
X = zeros(size(t_field,1), n);
Xdot = X;
dt = mean(diff(t_field));
for i = 1:n
    temporalAmplitude = V(:,i) * S(i,i);
    f = temporalAmplitude.';
    fdot = gradient(f,dt);
    X(:,i) = f';
    Xdot(:,i) = fdot';
end

polyorder = 2;
lambda = 0.0135;

Theta = poolData(X, n, polyorder);
Xi = sparsifyDynamics(Theta, Xdot, lambda, n);
 
fdot_sindy = Theta*Xi;
 
figure;
for j = 1:n
    f_sindy= cumtrapz(t_field, fdot_sindy(:,j));
    subplot(2,6,j);
    plot(t_field, fdot, 'LineWidth', 2); hold on
    plot(t_field, fdot_sindy(:,j), '--k', 'LineWidth',2);
    xlabel('Time(t)');
    ylabel('Derivative of f(t)');
    legend('True Derivative', 'SINDY Approximation');
    grid on
    title([num2str(j), 'th temporal amplitude']);

    subplot(2,6,j+6)
    plot(t_field, X(:,j), 'LineWidth', 2); hold on;
    plot(t_field, f_sindy + X(1,j), '--k', 'LineWidth', 2);
    xlabel('Time(t)');
    ylabel('f(t)')
    legend('True Function', 'SINDY derived');
    grid on
end

disp(Xi)
