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

%% Plotting diagonal elements to check the "energy content"
% Extract the diagonal elements (singular values) from S
singularValues = diag(S).^2; % Square of singular values

% Plot squared singular values on a linear-log scale
figure;
semilogy(singularValues, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'c');
grid on; % Add grid for better readability

% Add labels and title
xlabel('Mode index');
ylabel('Squared Singular Values');
xlim([0,10])
title('Energy Captured by Each Mode');

%% Plotting fractional energy content
% Calculate squared singular values
squaredSingularValues = diag(S).^2;

% Calculate the cumulative energy captured by each mode
cumulativeEnergy = cumsum(squaredSingularValues);

% Normalize by the total energy to get the fraction of the total energy
totalEnergy = sum(squaredSingularValues);
energyFraction = cumulativeEnergy / totalEnergy;

% Find the number of modes needed to capture 95% of the energy
threshold = 0.95;
modesFor95Energy = find(energyFraction >= threshold, 1, 'first');

% Plot the energy fraction of each mode
figure;
plot(energyFraction, 'o-', 'MarkerSize', 6, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
hold on;

% Plot a horizontal line at 95% energy capture
yline(threshold, '--', 'LineWidth', 2, 'Color', 'b');
hold off;

% Mark the point where 95% of the energy is captured
hold on;
plot(modesFor95Energy, energyFraction(modesFor95Energy), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
hold off;

% Add grid, labels, and title
grid on;
xlabel('Mode index');
ylabel('Cumulative Energy Fraction');
title('Cumulative Energy Content of Each Mode');
legend('Energy Fraction', '95% Energy', '95% Threshold Mode', 'Location', 'best');

% Adjust x-axis to focus on significant modes
xlim([1, modesFor95Energy + 10]);

%% Plotting first 6 spatial modes

% Number of modes to plot
numModes = 6;

% Adjust the layout
numRows = 3; % Increase the number of rows to spread out the plots
numCols = 4; % Adjust the number of columns as needed

% Create a large figure
figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.85]); % Use a large portion of the screen

% Total number of spatial points for each component
totalPoints = nx * ny;

chosen_colormap = 'viridis'; % Example: 'jet', 'parula', 'viridis'
contour_levels = 100; % Adjusted for smoother gradients
caxis_limits = [-0.1, 0.1]; % Adjusted range for enhanced contrast
%[~,hc] = contourf(x, y, uy_approx.', linspace(-1.5,1.5,80));
%    set(hc, 'LineStyle','none')
%    caxis([-0.1,0.1]); % Apply the chosen color axis limits
%    hold on;

for i = 1:numModes
    % Extract the mode for ux
    ux_mode = U(1:totalPoints, i); % First half for ux
    uy_mode = U(totalPoints+1:end, i); % Second half for uy
    
    % Reshape the modes to the original grid size
    ux_mode_reshaped = reshape(ux_mode, ny, nx);
    uy_mode_reshaped = reshape(uy_mode, ny, nx);
    
    % Plotting ux_mode
    subplot(numRows, numCols, i); % Adjust subplot positioning
    %[~,hc] = contourf(x,y,ux_mode_reshaped, linspace(-1.5,1.5,80));
    %set(hc, 'LineStyle','none');
    %caxis([-0.005,0.005])
    contourf(x,y,ux_mode_reshaped,20, 'LineStyle','none');
    colorbar;
    axis equal tight;
    title(['ux Mode ', num2str(i)]);
    
    % Ensure there's space for uy plots; adjust index calculation if necessary
    subplot(numRows, numCols, i+numModes); % Position uy plots in the layout
    %[~,hc] = contourf(x,y, uy_mode_reshaped, linspace(-1.5,1.5,80));
    %set(hc, 'LineStyle','none');
    %caxis([-0.1,0.1])
    contourf(x, y, uy_mode_reshaped, 20, 'LineStyle', 'none');
    colorbar;
    axis equal tight;
    title(['uy Mode ', num2str(i)]);
end

% Overall title for the subplot figure
sgtitle('First 6 Spatial Modes for ux and uy');

%% Spatial Mode Test

% Assuming U, totalPoints, nx, ny, dataMean are defined correctly
totalPoints = nx * ny

U_ux = U(1:totalPoints, :);
U_uy = U(totalPoints+1:end, :);

MM = 0.01;
v = -1:0.1:1;
v(11)=[];

figure;
for k = 1:6
    subplot(2, 3, k);
    contourf(x, y, transpose(reshape(U_ux(:,k), nx, ny)), MM * v);
    caxis([-MM MM]);
    colorbar;
end
sgtitle('ux spatial modes')

figure;
for k = 1:6
    subplot(2, 3, k);
    contourf(x, y, transpose(reshape(U_uy(:,k), nx, ny)), MM * v);
    caxis([-MM MM]);
    colorbar;
end
sgtitle('uy spatial modes')





%% Plotting first 6 temporal amplitudes

% Number of modes to plot
numModes = 6

% Create a new figure
figure;
hold on;

% Loop through the first 6 modes
for i = 1:numModes
    % Calculate the temporal amplitude for the i-th mode
    % V(:, i) is the temporal mode, and S(i, i) is the corresponding singular value
    temporalAmplitude = V(:, i) * S(i, i);
    
    % Plot the temporal amplitude
    plot(t_field, temporalAmplitude, 'LineWidth', 2, 'DisplayName', ['Mode ', num2str(i)]);
end

% Improve plot appearance
xlabel('Time');
ylabel('Temporal Amplitude');
title('First 6 Temporal Amplitudes');
legend show;
grid on;
hold off;


%% Making a video
% Setup for video creation
videoFilename = 'FluidFlowReconstruction.avi';
v = VideoWriter(videoFilename, 'Uncompressed AVI'); % Choose your desired format
v.FrameRate = 10; % Adjust frame rate as needed
open(v);

% Loop through each time step
for t = 1:nt
    % Create a new figure for each frame to ensure consistent sizing
    fig = figure('Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Visible', 'off');
    
    % Original flow plot
    subplot(1, 5, 1);
    imagesc(x, y, squeeze(uy(:,:,t))');
    axis equal tight;
    title('Original Flow');
    colorbar;
    
    % Ranks for reconstruction
    ranks = [2, 4, 10, 20];
    
    for i = 1:length(ranks)
        r = ranks(i);
        % Reconstruct flow using the first 'r' modes
        data_approx = U(:,1:r) * S(1:r,1:r) * V(:,1:r)' + repmat(dataMean, 1, nt); % Add back the mean
        uy_approx = reshape(data_approx(totalPoints+1:end, t), ny, nx);
        
        % Plot reconstruction
        subplot(1, 5, i+1);
        imagesc(x, y, uy_approx');
        axis equal tight;
        title(['Reconstruction r=', num2str(r)]);
        colorbar;
    end
    
    drawnow;
    
    % Capture and write frame
    frame = getframe(fig);
    writeVideo(v, frame);
    
    % Close the figure to free up memory
    close(fig);
end

close(v);


%% Video

% Choose the rank for reconstruction
r = 401;

% Prepare the video file
videoFilename = 'POD_Reconstruction_r401.avi';
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
    data_approx = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)' + repmat(dataMean, 1, 1);
    uy_approx = reshape(data_approx(nx*ny+1:end, t), nx, ny);
    
    % Clear the previous plot
    clf;
   
    % Generate the plot for the current time step
    [~,hc] = contourf(x, y, uy_approx.', linspace(-1.5,1.5,80));
    set(hc, 'LineStyle','none')
    caxis([-0.1,0.1]); % Apply the chosen color axis limits
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

%% Original Video

% Prepare the video file for the original dataset
videoFilenameOriginal = 'Original_Dataset.avi';
vOriginal = VideoWriter(videoFilenameOriginal);
vOriginal.FrameRate = 60; % Match the frame rate of the reconstruction video for consistency
open(vOriginal);

% Set up the figure for video frames
figure;

chosen_colormap = 'viridis'; % Keep consistent visualization settings
contour_levels = 100; % Adjusted for smoother gradients
caxis_limits = [-0.1, 0.1]; % Adjusted range for enhanced contrast

% Loop over each time step to visualize the original data
for t = 1:nt
    % Access the original data for the current timestep
    % Assuming 'uy_original' is your original uy component data
    uy_original = reshape(uy(:, :, t), nx, ny); % Reshape if necessary
    
    % Clear the previous plot
    clf;
    
    % Generate the plot for the current time step using the original data
    [~,hc] = contourf(x, y, uy_original.', linspace(-1.5, 1.5, 80));
    set(hc, 'LineStyle', 'none');
    caxis(caxis_limits); % Apply the chosen color axis limits
    hold on;
    % Plot airfoil contour if applicable
    plot(xa(:,t), ya(:,t), 'k-', 'LineWidth', 1.5);
    xlim([min(x) - 0.1, max(x) + 0.1]);
    ylim([-1, 1]);
    title(sprintf('Original Fluid Flow, Time Step = %d', t), 'Interpreter', 'Latex');
    xlabel('$x/c$', 'Interpreter', 'Latex');
    ylabel('$y/c$', 'Interpreter', 'Latex');
    colorbar;
    axis equal;
    axis tight;
    
    % Capture the frame for the video
    frame = getframe(gcf);
    writeVideo(vOriginal, frame);
end

% Close the video file
close(vOriginal);

% Close the figure
close(gcf);

%% Subplots

% Video setup
videoFilename = 'Flow_Reconstruction_Comparison_MultiRow.avi';
v = VideoWriter(videoFilename, 'Motion JPEG AVI'); % Using MJPEG for compatibility and quality
v.FrameRate = 60; % High frame rate for smoother playback
v.Quality = 100; % Maximum quality
open(v);

% Define the ranks for reconstruction
ranks = [2, 4, 6, 8, 10];
totalPlots = length(ranks) + 1; % Total number of plots including the original

% Determine the layout for subplots
numRows = 2; % Number of rows to display the plots
numCols = ceil(totalPlots / numRows); % Calculate columns needed

% Set up the figure for video frames, specifying high resolution
figure('Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Visible', 'off');

% Visualization settings
chosen_colormap = 'viridis';
contour_levels = 100;
caxis_limits = [-0.1, 0.1];

% Loop over each time step for the reconstruction
for t = 1:nt
    clf; % Clear figure for new frame
    
    % Original flow plot
    subplot(numRows, numCols, 1); % Position original flow in the first subplot

    uy_original = reshape(uy(:, :, t), nx, ny); % Reshape if necessary
    
    % Generate the plot for the current time step using the original data
    [~,hc] = contourf(x, y, uy_original.', linspace(-1.5, 1.5, 80));
    set(hc, 'LineStyle', 'none');
    caxis(caxis_limits); % Apply the chosen color axis limits
    hold on;
    % Plot airfoil contour if applicable
    plot(xa(:,t), ya(:,t), 'k-', 'LineWidth', 1.5);
    xlim([min(x) - 0.1, max(x) + 0.1]);
    ylim([-1, 1]);
    title(sprintf('Original Fluid Flow, Time Step = %d', t), 'Interpreter', 'Latex');
    xlabel('$x/c$', 'Interpreter', 'Latex');
    ylabel('$y/c$', 'Interpreter', 'Latex');
    colorbar;
    axis equal;
    axis tight;
    
    % Loop over each rank for reconstructions
    for i = 1:length(ranks)
        r = ranks(i);
        data_approx = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)' + repmat(dataMean, 1, nt);
        uy_approx = reshape(data_approx(nx*ny+1:end, t), nx, ny);
        
        % Calculate subplot index for current reconstruction
        subplotIndex = i + 1; % Adjust subplot index based on the current rank
        subplot(numRows, numCols, subplotIndex);
        

        % Generate the plot for the current time step
        [~,hc] = contourf(x, y, uy_approx.', linspace(-1.5,1.5,80));
        set(hc, 'LineStyle','none')
        caxis([-0.1,0.1]); % Apply the chosen color axis limits
        colorbar;
        title(['Reconstruction r=', num2str(r)]);
        axis equal tight;
    end
    
    drawnow; % Update figure
    
    % Capture the frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);

% Optionally, close the figure
close(gcf);
