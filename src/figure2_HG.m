% parameter definition
omega0 = 1; % beam waist radius
zR = 10;    % Rayleigh length
lambda = 0.6328e-6; % wavelength (meters)
f = 0.1;    % lens focal length (meters)
k = 2 * pi / lambda; % wave number

% define functions
omega = @(z) omega0 * sqrt(1 + (z / zR)^2); % beam width
gouy = @(z) atan(z / zR);                   % Gouy phase

% HG mode function
HG_mode = @(n, m, x, y, z) ...
    sqrt(2^(-(n + m)) / (pi * factorial(n) * factorial(m))) .* ...
    hermiteH(n, sqrt(2) * x ./ omega(z)) .* ...
    hermiteH(m, sqrt(2) * y ./ omega(z)) .* ...
    exp(-(x.^2 + y.^2) ./ omega(z).^2) .* ...
    exp(1i * (n + m + 1) * gouy(z));

% define coordinate grid
x = linspace(-5, 5, 300); % coordinate range
y = linspace(-5, 5, 300);
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);
z = 10; % specify z position

% define figure size
figure('Position', [100, 100, 1200, 900]);

% parameter settings (select specific modes)
n_values = [2, 3, 2]; % HG mode x order n
m_values = [0, 1, 3]; % HG mode y order m

% -----------------
% HG mode plotting
% -----------------
% First row: HG Amplitude
for i = 1:length(n_values)
    n = n_values(i);
    m = m_values(i);
    amplitude = abs(HG_mode(n, m, X, Y, z)); % amplitude
    subplot(4, length(n_values), i + 3 * length(p_values)); % first row
    imagesc(x, y, amplitude);
    axis off; axis square;
    title(['HG_{', num2str(n), num2str(m), '} (Amplitude)']);
    colormap hot;
end

% Second row: HG Phase
for i = 1:length(n_values)
    n = n_values(i);
    m = m_values(i);
    phase = angle(HG_mode(n, m, X, Y, z)); % phase
    subplot(4, length(n_values), i + 4 * length(p_values)); % second row
    imagesc(x, y, phase);
    axis off; axis square;
    title(['HG_{', num2str(n), num2str(m), '} (Phase)']);
    colormap hsv;
end

% Third row: HG Phase-only hologram
for i = 1:length(n_values)
    n = n_values(i);
    m = m_values(i);
    % HG mode phase
    mode_phase = angle(HG_mode(n, m, X, Y, z));
    % add grating phase
    grating_phase = 2 * pi * (X / x0 + Y / y0);
    % phase hologram
    phase_hologram = mod(mode_phase + grating_phase, 2 * pi);
    subplot(4, length(n_values), i + 5 * length(p_values)); % third row
    imagesc(x, y, phase_hologram);
    axis off; axis square;
    title(['HG_{', num2str(n), num2str(m), '} (Phase Hologram)']);
    colormap gray;
end
