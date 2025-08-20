% parameter settings
omega0 = 1; % beam waist radius
zR = 10;    % Rayleigh length
lambda = 0.532e-6; % wavelength (meters)
f = 0.25e-6;    % lens focal length (meters)
k = 2 * pi / lambda; % wave number

% define functions
omega = @(z) omega0 * sqrt(1 + (z / zR)^2); % beam width
gouy = @(z) atan(z / zR);                   % Gouy phase

% LG mode function
LG_mode = @(p, l, rho, phi, z) ...
    sqrt(2 * factorial(p) / (pi * factorial(abs(l) + p))) .* ...
    (sqrt(2) * rho ./ omega(z)).^abs(l) .* ...
    laguerreL(p, abs(l), 2 * rho.^2 ./ omega(z).^2) .* ...
    exp(-rho.^2 ./ omega(z).^2) .* ...
    exp(1i * (2 * p + abs(l) + 1) * gouy(z)) .* ...
    exp(-1i * l * phi);

% define coordinate grid
x = linspace(-5, 5, 300); % coordinate range
y = linspace(-5, 5, 300);
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
z = 10; % specify z position

% define figure size
figure('Position', [100, 100, 1200, 600]);

% parameter settings (select specific modes)
p_values = [2, 2, 1]; % LG mode radial order p
l_values = [0, 1, 2]; % LG mode azimuthal order l

% % first row: Amplitude
% for i = 1:length(p_values)
%     p = p_values(i);
%     l = l_values(i);
%     amplitude = abs(LG_mode(p, l, rho, phi, z)); % amplitude
%     subplot(3, length(p_values), i); % first row
%     imagesc(x, y, amplitude);
%     axis off; axis square;
%     title(['LG^{', num2str(l), '}_{', num2str(p), '} (Amplitude)']);
%     colormap parula;
% end
%
% % second row: Phase
% for i = 1:length(p_values)
%     p = p_values(i);
%     l = l_values(i);
%     phase = angle(LG_mode(p, l, rho, phi, z)); % phase
%     subplot(3, length(p_values), i + length(p_values)); % second row
%     imagesc(x, y, phase);
%     axis off; axis square;
%     title(['LG^{', num2str(l), '}_{', num2str(p), '} (Phase)']);
%     colormap parula;
% end

% third row: Phase-only hologram
x0 = 50; y0 = 50; % grating period parameters
for i = 1:length(p_values)
    p = p_values(i);
    l = l_values(i);
    % LG mode phase
    mode_phase = angle(LG_mode(p, l, rho, phi, z));
    % % add grating phase
    % grating_phase = 2 * pi * (X / x0 + Y / y0);
    % phase hologram
    phase_hologram = mod(mode_phase, 2 * pi);
    subplot(3, length(p_values), i + 2 * length(p_values)); % third row
    imagesc(x, y, phase_hologram);
    axis off; axis square;
    title(['LG^{', num2str(l), '}_{', num2str(p), '} (Phase Hologram)']);
    colormap gray;
end
