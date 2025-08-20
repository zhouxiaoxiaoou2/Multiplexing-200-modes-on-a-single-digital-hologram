clc; clear; close all;
rng(1); % set random seed for reproducibility

% ------ Define Laguerre-Gaussian mode function ------
omega0 = 1; zR = 5; % beam waist and Rayleigh length
omega = @(z) omega0 * sqrt(1 + (z / zR)^2);
LG_mode = @(p, l, rho, phi, z) ...
    sqrt(2 * factorial(p) / (pi * factorial(abs(l) + p))) .* ...
    (sqrt(2) * rho ./ omega(z)).^abs(l) .* ...
    laguerreL(p, abs(l), 2 * rho.^2 ./ omega(z).^2) .* ...
    exp(-rho.^2 ./ omega(z).^2) .* ...
    exp(1i * (2 * p + abs(l) + 1) .* atan(z / zR));

% ------ Define coordinate system ------
x = linspace(-5, 5, 100); % spatial coordinates
y = linspace(-5, 5, 100);
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
z = 0; % beam focus position

% ------ Set different multiplexing mode numbers ------
num_modes_list = [3, 50, 75, 100]; 
correlation_phase_list = zeros(size(num_modes_list));
correlation_cam_list = zeros(size(num_modes_list));

% ------ Theoretical mode (fixed LG_2^1 mode) ------
p_theory = 2; l_theory = 1; 
theoretical_mode = abs(LG_mode(p_theory, l_theory, rho, phi, z)); % Theoretical mode amplitude

% ------ Progress bar ------
h = waitbar(0, 'Processing...');

for idx = 1:length(num_modes_list)
    num_modes = num_modes_list(idx); % current multiplexing mode number
    waitbar(idx / length(num_modes_list), h, sprintf('Processing %d modes...', num_modes));

    % ------ Generate multiple random modes ------
    mixed_modes = zeros(size(X));
    for m = 1:num_modes
        p = randi([0, 4]); % randomly select p
        l = randi([-4, 4]); % randomly select l
        mixed_modes = mixed_modes + abs(LG_mode(p, l, rho, phi, z));
    end
    mixed_modes = mixed_modes / max(mixed_modes(:)); % normalize

    % ------ Phase-only modulation ------
    phase_noise = exp(1i * 0.1 * randn(size(mixed_modes))); % control phase noise
    LG_phase_only = abs(ifft2(fft2(mixed_modes) .* phase_noise)).^2; % inverse transform to get field intensity

    % ------ CAM modulation ------
    amp_mod = 1 + 0.05 * randn(size(mixed_modes)); % control amplitude noise
    LG_CAM = abs(ifft2(fft2(mixed_modes) .* amp_mod)).^2;

    % ------ Compute normalized correlation ------
    correlation_phase_list(idx) = corr2(LG_phase_only / max(LG_phase_only(:)), theoretical_mode / max(theoretical_mode(:)));
    correlation_cam_list(idx) = corr2(LG_CAM / max(LG_CAM(:)), theoretical_mode / max(theoretical_mode(:)));
    %
    % ------ Display results ------
    figure;
    subplot(2,3,1); imagesc(theoretical_mode); axis off; colormap('parula'); title('Theory');
    subplot(2,3,2); imagesc(angle(phase_noise)); axis off; colormap('parula'); title(['Phase-only Noise']);
    subplot(2,3,3); imagesc(LG_phase_only); axis off; colormap('parula'); title(['Recovered Phase-only, C=', num2str(correlation_phase_list(idx))]);
    subplot(2,3,4); imagesc(theoretical_mode); axis off; colormap('parula'); title('Theory');
    subplot(2,3,5); imagesc(amp_mod); axis off; colormap('parula'); title(['CAM Noise']);
    subplot(2,3,6); imagesc(LG_CAM); axis off; colormap('parula'); title(['Recovered CAM, C=', num2str(correlation_cam_list(idx))]);
end

close(h); % close progress bar

% ------ Plot Phase-only correlation ------
figure;
plot(num_modes_list, correlation_phase_list, 'ro-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
yline(0.8, '--b', 'Threshold', 'LineWidth', 1.5);
ylim([0.3 1.2]); % set y-axis limits
xlabel('Number of Modes');
ylabel('Correlation');
title('Phase-only');
legend('Phase-only', 'Threshold');
grid on;

% ------ Plot CAM correlation ------
figure;
plot(num_modes_list, correlation_cam_list, 'ro-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
yline(0.8, '--b', 'Threshold', 'LineWidth', 1.5);
ylim([0.3 1.2]); % set y-axis limits
xlabel('Number of Modes');
ylabel('Correlation');
title('Complex Amplitude Modulation');
legend('CAM', 'Threshold');
grid on;