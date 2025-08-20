% set figure properties
figure;
set(gcf, 'Position', [100, 100, 1800, 300]);

% define grid
x = linspace(-5, 5, 300); % use appropriate resolution
y = linspace(-5, 5, 300);
[X, Y] = meshgrid(x, y);
[theta, r] = cart2pol(X, Y);

% define carrier frequencies for generating tilted fringes
frequency1 = 20;  % mode 1 fringe frequency
frequency2 = 15;  % mode 2 fringe frequency
frequency3 = 15;  % mode 3 fringe frequency

angle_deg1 = 225;  % mode 1 fringe angle
angle_deg2 = 180;  % mode 2 fringe angle
angle_deg3 = 270; % mode 3 fringe angle

% define three carriers
angle_rad1 = deg2rad(angle_deg1);
carrier1 = exp(1i * frequency1 * (cos(angle_rad1) * X + sin(angle_rad1) * Y));

angle_rad2 = deg2rad(angle_deg2);
carrier2 = exp(1i * frequency2 * (cos(angle_rad2) * X + sin(angle_rad2) * Y));

angle_rad3 = deg2rad(angle_deg3);
carrier3 = exp(1i * frequency3 * (cos(angle_rad3) * X + sin(angle_rad3) * Y));

% define LG mode amplitude and phase distribution functions
LG_amplitude = @(p, l, r) abs((r.^l) .* exp(-r.^2) .* laguerreL(p, l, r.^2));
LG_phase = @(p, l, theta, r) angle((r.^l) .* exp(-r.^2) .* laguerreL(p, l, r.^2) .* exp(1i * l * theta));

% define CAM hologram functions
LG_CAM_hologram1 = @(p, l, theta, r) LG_amplitude(p, l, r) .* exp(1i * LG_phase(p, l, theta, r)) .* carrier1;
LG_CAM_hologram2 = @(p, l, theta, r) LG_amplitude(p, l, r) .* exp(1i * LG_phase(p, l, theta, r)) .* carrier2;
LG_CAM_hologram3 = @(p, l, theta, r) LG_amplitude(p, l, r) .* exp(1i * LG_phase(p, l, theta, r)) .* carrier3;
%
% generate LG mode CAM holograms
subplot(1, 6, 1);
imagesc(real(LG_CAM_hologram1(1, 1, theta, r))); % LG_1^1模式
axis square;
axis off;
title('LG_1^1 CAM Hologram');

subplot(1, 6, 2);
imagesc(real(LG_CAM_hologram2(2, 0, theta, r))); % LG_0^2模式
axis square;
axis off;
title('LG_2^0 CAM Hologram');

subplot(1, 6, 3);
imagesc(real(LG_CAM_hologram3(0, 2, theta, r))); % LG_2^0模式
axis square;
axis off;
title('LG_0^2 CAM Hologram');
% 
multiplexed_hologram = LG_CAM_hologram1(1, 1, theta, r) + LG_CAM_hologram2(2, 0, theta, r) + LG_CAM_hologram3(0, 2, theta, r);
% 
subplot(1, 6, 4);
imagesc(real(multiplexed_hologram)); 
axis square;
axis off;
title('Multiplexed Hologram');

% adjust colormap
colormap(gray); % use gray colormap


%%%Fourier Transform

% simulate Fourier lens
    Fourier = fftshift(fft2(ifftshift(multiplexed_hologram)));

    % define spatial filter
    [fx, fy] = meshgrid(linspace(-1, 1, size(Fourier, 1)));
    filter_radius = 0.1;
    spatial_filter = sqrt(fx.^2 + fy.^2) < filter_radius;

    % apply spatial filter and reconstruct
    filtered = Fourier .* spatial_filter;
    amplitude = abs(ifftshift(ifft2(fftshift(filtered))));


% Step 2: Display the Fourier Transform of the multiplexed hologram (amplitude only)
figure;
imagesc(amplitude); % display amplitude
colormap parula;
colorbar;
title('Fourier Transform of Multiplexed Hologram (Amplitude Only)');
axis image;

