clc; clear all; 

%Parameters for the transfer function
zeta1 = 0.03;
zeta2 = 0.03;
zeta3 = 0.042;
zeta4 = 0.025;
zeta5 = 0.032;
k = 5;
tau = 10^-4;
w1 = 2*pi*2.4*1000;
w2 = 2*pi*2.6*1000;
w3 = 2*pi*6.5*1000;
w4 = 2*pi*8.3*1000;
w5 = 2*pi*9.3*1000;
nin = 1; nout = 1;

%Transfer function
s = tf('s');
num = k*w2^2*w3^2*w5^2*(s^2 + 2*zeta1*w1*s + w1^2)*(s^2 + 2*zeta4*w4*s + w4^2)*exp(-s*tau);
denom = w1^2*w4^2*(s^2 + 2*zeta2*w2*s + w2^2)*(s^2 + 2*zeta3*w3*s + w3^2)*(s^2 + 2*zeta5*w5*s + w5^2);
H = num/denom;

%Getting impulse response
[y, t] = impulse(H);

% Eigensystem Realization Algorithm
H0 = hankel(y(1:end-1));
H1 = hankel(y(2:end));

[U, Sigma, V] = svd(H0, 'econ'); %SVD
r = 15;%<< order of the desired reduced system
Ur = U(:,1:r); %<< truncating SVD components
Sigma_r = Sigma(1:r, 1:r);
Vr = V(:, 1:r);

Ip = zeros(size(Vr,1),1); Ip(1) = 1;
Iq = zeros(1,size(Ur,1)); Iq(1) = 1;
A = Sigma_r^(-0.5) * Ur' * H1 * Vr * Sigma_r^(-0.5);
B = Sigma_r^(0.5) * Vr' * Ip;
C = Iq * Ur * Sigma_r^(0.5);
D = y(1);

%Frequency response comparison
figure;
bode(H);
grid on;
hold on;

dt = t(2)-t(1);

opts = bodeoptions('cstprefs');
opts.PhaseWrapping = 'on';
bode(dt*ss(A,B,C,0,ts = dt),opts)

% Testing the robustness of ERA
noisy_y = y + 0.8 * randn(size(y));
H0 = hankel(noisy_y(1:end-1));
H1 = hankel(noisy_y(2:end));

[U, Sigma, V] = svd(H0, 'econ'); %SVD
r = 15;%<< order of the desired reduced system
Ur = U(:,1:r); %<< truncating SVD components
Sigma_r = Sigma(1:r, 1:r);
Vr = V(:, 1:r);

Ip = zeros(size(Vr,1),1); Ip(1) = 1;
Iq = zeros(1,size(Ur,1)); Iq(1) = 1;
A = Sigma_r^(-0.5) * Ur' * H1 * Vr * Sigma_r^(-0.5);
B = Sigma_r^(0.5) * Vr' * Ip;
C = Iq * Ur * Sigma_r^(0.5);
D = y(1);

opts = bodeoptions('cstprefs');
opts.PhaseWrapping = 'on';
bode(dt*ss(A,B,C,0,ts = dt),opts)
legend('original','ERA', 'ERA with noise')
