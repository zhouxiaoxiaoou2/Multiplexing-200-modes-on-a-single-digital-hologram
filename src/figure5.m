% parameter settings
omega0 = 1.; % beam waist radius

% define LG mode function
LG_mode = @(p, l, rho, phi) ...
    (1 / omega0) * sqrt((2 * factorial(p)) / (pi * factorial(l + p))) * ...
    (sqrt(2) * rho / omega0).^abs(l) .* laguerreL(p, abs(l), (2 * rho.^2) / omega0^2) .* ...
    exp(-rho.^2 / omega0^2) .* exp(-1i * l * phi);

% define HG mode function
HG_mode = @(n, m, x, y) ...
    (1 / omega0) * sqrt((2^(-n - m) * factorial(n) * factorial(m)) / pi) * ...
    hermiteH(n, sqrt(2) * x / omega0) .* hermiteH(m, sqrt(2) * y / omega0) .* ...
    exp(-(x.^2 + y.^2) / omega0^2);

% define coordinate grid
x = linspace(-5, 5, 300);
y = linspace(-5, 5, 300);
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2);
phi = atan2(Y, X);

% define carrier frequency and angle
frequency = 10; % fringe frequency
angle_deg = 45; % fringe angle
angle_rad = deg2rad(angle_deg);
carrier = exp(1i * frequency * (cos(angle_rad) * X + sin(angle_rad) * Y));

% mode information
modes = {'LG_3^2', 'LG_4^4', 'HG_{44}', 'HG_{73}'};
C_values = zeros(1, length(modes)); % for storing correlation coefficients

% custom correlation computation function
function C = compute_correlation(A, B)
    A_mean = mean(A(:));
    B_mean = mean(B(:));
    numerator = sum(sum((A - A_mean) .* (B - B_mean)));
    denominator = sqrt(sum(sum((A - A_mean).^2)) * sum(sum((B - B_mean).^2)));
    C = numerator / denominator;
end

% iterate over all modes
for i = 1:length(modes)
    switch modes{i}
        case 'LG_3^2'
            theory = abs(LG_mode(3, 2, rho, phi)).^2;
            hologram = abs(LG_mode(3, 2, rho, phi) .* carrier).^2;
        case 'LG_4^4'
            theory = abs(LG_mode(4, 4, rho, phi)).^2;
            hologram = abs(LG_mode(4, 4, rho, phi) .* carrier).^2;
        case 'HG_{44}'
            theory = abs(HG_mode(4, 4, X, Y)).^2;
            hologram = abs(HG_mode(4, 4, X, Y) .* carrier).^2;
        case 'HG_{73}'
            theory = abs(HG_mode(7, 3, X, Y)).^2;
            hologram = abs(HG_mode(7, 3, X, Y) .* carrier).^2;
    end

    % simulate Fourier lens
    Fourier = fftshift(fft2(ifftshift(hologram)));

    % define spatial filter
    [fx, fy] = meshgrid(linspace(-1, 1, size(Fourier, 1)));
    filter_radius = 0.1;
    spatial_filter = sqrt(fx.^2 + fy.^2) < filter_radius;

    % apply spatial filter and reconstruct
    filtered = Fourier .* spatial_filter;
    amplitude = abs(ifftshift(ifft2(fftshift(filtered))));

    % compute correlation coefficient
    C_values(i) = compute_correlation(theory, amplitude);

    % plot mode distribution
    figure;
    subplot(2, 2, 1);
    imagesc(theory); colorbar; title([modes{i}, ' Theoretical Intensity']);
    subplot(2, 2, 2);
    imagesc(amplitude); colorbar; title([modes{i}, ' Experimental Intensity']);

    % cross-section plot
    middle_row = round(size(theory, 1) / 2);
    subplot(2, 2, [3, 4]);
    plot(x, theory(middle_row, :), 'g-', 'DisplayName', 'Theory'); hold on;
    plot(x, amplitude(middle_row, :), 'b--', 'DisplayName', 'Exp');
    legend;
    title(['Cross-section Analysis - ', modes{i}]);
    xlabel('Position');
    ylabel('Intensity');
end

% print correlation coefficients
for i = 1:length(modes)
    fprintf('%s Correlation Coefficient C = %.2f\n', modes{i}, C_values(i));
end
