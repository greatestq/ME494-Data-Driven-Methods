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

BaseAngle = 30; % options are 25 and 30 deg

Freqs = [0.05,0.25,0.5]; % choose pitching frequencies to plot
% must be from the set 0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5

FreqLabels = {'0p05','0p25','0p5'}; % for loading files, must match selected Freqs

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

%% Truncation and Plotting
figure
r = [2 4 6 8 10 20 50];
s = size(data,1)/2 + 1;
%s = nx*ny + 1;

%Contour Parameters
chosen_colormap = 'viridis'; %Example: 'jet', 'parula', 'viridis'
contour_levels = 100; %Increase for smoother gradients
caxis_limits = [-0.1, 0.1]; %Narrower range for enhanced contrast

for i = 1:length(r)
    subplot(length(r),1,i)
    data_approx = U(:,1:r(i))*S(1:r(i),1:r(i))*V(:,1:r(i))';
    uy_approx = reshape(data_approx(s:end, tstep), nx, ny);
    %[~,hc] = contourf(x,y,squeeze(uy(:,:,tstep)).',linspace(-1.5,1.5,80));
    [~,hc] = contourf(x,y,uy_approx.', linspace(-1.5,1.5,80));
    set(hc,'LineStyle','none');
    caxis([-0.5,0.5]); % plot with same range for each case
    hold on
    % plot airfoil coordinates at appropriate timestep
    plot(xa(:,tstep),ya(:,tstep),'k-')
end

%% POD analysis based on "short" simulation
% Number of time steps to include in the analysis
nt_short = 100;

%Timestep of interest
tstep = 50;

% Restrict the data to the first 100 time steps
ux_short = ux(:, :, 1:nt_short);
uy_short = uy(:, :, 1:nt_short);

% Reshape the restricted data for SVD
uxreshape_short = reshape(ux_short, nx*ny, nt_short);
uyreshape_short = reshape(uy_short, nx*ny, nt_short);
data_short = [uxreshape_short; uyreshape_short];

% Subtract the mean if meanSub is enabled
if meanSub
    dataMean_short = mean(data_short, 2);
    data_short = data_short - dataMean_short*ones(1, nt_short);
end

% Perform SVD on the restricted data
[U_short, S_short, V_short] = svd(data_short, 'econ');

% Define the number of modes you want to use for reconstruction
r = [2, 4, 6, 8, 10, 20, 50]; % Example set of ranks for reconstruction

% Adjusted contour parameters for clearer visualization
chosen_colormap = 'viridis'; %Example: 'jet', 'parula', 'viridis'
contour_levels = 100; %Increase for smoother gradients
caxis_limits = [-0.1, 0.1]; %Narrower range for enhanced contrast

figure;
s = nx*ny + 1; % Starting index for uy data in the reshaped matrix

for i = 1:length(r)
    subplot(length(r), 1, i);
    % Reconstruct the flow field using the first r(i) modes
    data_approx_short = U_short(:, 1:r(i)) * S_short(1:r(i), 1:r(i)) * V_short(:, 1:r(i))';
    % Separate and reshape the uy component of the reconstructed flow field
    uy_approx_short = reshape(data_approx_short(s:end, tstep), nx, ny);
    
    % Create the contour plot for the reconstructed uy component
    [~,hc] = contourf(x,y,uy_approx.', linspace(-1.5,1.5,80));
    set(hc,'LineStyle','none');
    caxis([-0.5,0.5]); % plot with same range for each case
    
    hold on;
    % Plot airfoil coordinates at the chosen timestep
    plot(xa(:,tstep), ya(:,tstep), 'k-', 'LineWidth', 1.5);
    
    title(sprintf('Reconstructed $u_y$ using %d modes (First 100 steps)', r(i)), 'Interpreter', 'Latex');
    xlabel('$x/c$', 'Interpreter', 'Latex');
    ylabel('$y/c$', 'Interpreter', 'Latex');
    
    if i == length(r) % Add a colorbar to the last subplot for reference
        colorbar;
    end
end

% Adjust overall figure properties if necessary
sgtitle('POD Reconstruction with First 100 Time Steps', 'Interpreter', 'Latex');


%% Energy Plot
figure
subplot(1,2,1), semilogy(diag(S),'k')
subplot(1,2,2), plot(cumsum(diag(S))/sum(diag(S)),'k')

%% Video

% Choose the rank for reconstruction
r = 10;

% Prepare the video file
videoFilename = 'POD_Reconstruction_r10.avi';
v = VideoWriter(videoFilename);
v.FrameRate = 60; % Increased frame rate to speed up the video
open(v);

% Set up the figure for video frames
figure;

chosen_colormap = 'viridis'; % Example: 'jet', 'parula', 'viridis'
contour_levels = 100; % Adjusted for smoother gradients
caxis_limits = [-0.1, 0.1]; % Adjusted range for enhanced contrast

% Loop over each time step for the reconstruction
for t = 1:nt
    % Reconstruct the flow field at time step t
    data_approx = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';
    uy_approx = reshape(data_approx(nx*ny+1:end, t), nx, ny);
    
    % Clear the previous plot
    clf;
   
    % Generate the plot for the current time step
    [~,hc] = contourf(x, y, uy_approx.', linspace(-1.5,1.5,80));
    set(hc, 'LineStyle','none')
    caxis([-0.5,0.5]); % Apply the chosen color axis limits
    hold on;
    % Ensure xa and ya are updated for each frame if they change over time
    % This assumes xa and ya are matrices with dimensions [numPoints x nt] where nt is the number of timesteps
    plot(xa(:,t), ya(:,t), 'k-', 'LineWidth', 1.5); % Animated airfoil contour for each timestep
    xlim([min(xa(:, t)) - 0.1, max(xa(:, t)) + 0.1]);
    ylim([-1, 1]);
    title(sprintf('Flow Field Reconstruction at r = %d, Time Step = %d', r, t), 'Interpreter', 'Latex');
    xlabel('$x/c$', 'Interpreter', 'Latex');
    ylabel('$y/c$', 'Interpreter', 'Latex');
    colorbar;
    axis equal;
    axis tight;
    
    % Capture the frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);

% Close the figure
close(gcf);

%% Pitching Motion test
% Prepare the figure
figure;

% Set the axis limits based on the entire dataset to keep a consistent view
x_limits = [min(xa(:)) - 0.001, max(xa(:)) + 0.001];
y_limits = [min(ya(:)) - 0.001, max(ya(:)) + 0.001];

% Loop through each time step
for t = 1:nt
    plot(xa(:,t), ya(:,t), 'k-', 'LineWidth', 1.5); % Plot the airfoil contour at timestep t
    axis equal; % Keep aspect ratio equal to ensure correct visualization of the airfoil shape
    xlim(x_limits); % Apply consistent x-axis limits
    ylim(y_limits); % Apply consistent y-axis limits
    title(sprintf('Airfoil Position at Time Step: %d', t));
    xlabel('x');
    ylabel('y');
    
    pause(0.1); % Pause to visually inspect each frame. Adjust as needed for faster or slower playback.
    
    % Optionally, clear the figure for the next timestep (uncomment the next line if needed)
    % clf;
end

