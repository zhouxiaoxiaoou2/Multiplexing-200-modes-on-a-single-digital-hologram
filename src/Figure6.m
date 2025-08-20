% parameter settings
omega0 = 1.0; % beam waist radius
x = linspace(-5, 5, 300); % define x coordinates
y = linspace(-5, 5, 300); % define y coordinates
[X, Y] = meshgrid(x, y); % create coordinate grid
rho = sqrt(X.^2 + Y.^2); % polar coordinate rho
phi = atan2(Y, X); % polar coordinate phi

% define LG mode function
LG_mode = @(p, l, rho, phi) ...
    (1 / omega0) * sqrt((2 * factorial(p)) / (pi * factorial(l + p))) * ...
    (sqrt(2) * rho / omega0).^abs(l) .* ...
    laguerreL(p, abs(l), (2 * rho.^2) / omega0^2) .* ...
    exp(-rho.^2 / omega0^2) .* exp(-1i * l * phi);

% define HG mode function
HG_mode = @(n, m, x, y) ...
    (1 / omega0) * sqrt((2^(-n - m) * factorial(n) * factorial(m)) / pi) .* ...
    hermiteH(n, sqrt(2) * x / omega0) .* ...
    hermiteH(m, sqrt(2) * y / omega0) .* ...
    exp(-(x.^2 + y.^2) / omega0^2);

% mode parameters
LG_params = {[3, 1], [0, 2], [2, 3], [4, 4]}; % LG mode parameters (p, l)
HG_params = {[1, 2], [3, 1], [2, 3], [4, 4]}; % HG mode parameters (n, m)

modes = {'LG_{3}^{1}', 'LG_{0}^{2}', 'LG_{2}^{3}', 'LG_{4}^{4}', ...
         'HG_{1,2}', 'HG_{3,1}', 'HG_{2,3}', 'HG_{4,4}'};

% initialize storage for images
theory_profiles = {};
phase_only_profiles = {};
cam_profiles = {};

% iterate over all modes
for k = 1:4
    % LG mode
    p = LG_params{k}(1);
    l = LG_params{k}(2);
    LG_theory = abs(LG_mode(p, l, rho, phi)).^2;
    
    % Phase-only Hologram
    phase_noise = exp(1i * rand(size(LG_theory))); % random phase noise
    LG_phase_only = abs(ifft2(fft2(LG_theory) .* phase_noise)).^2;

    % CAM Hologram
    amp_mod = 1 + 0.1 * rand(size(HG_theory)); % amplitude modulation
    LG_CAM = abs(ifft2(fft2(LG_theory) .* amp_mod)).^2;

    % save results
    theory_profiles{k} = LG_theory / max(LG_theory(:));
    phase_only_profiles{k} = LG_phase_only / max(LG_phase_only(:));
    cam_profiles{k} = LG_CAM / max(LG_CAM(:));
end

for k = 1:4
    % HG mode
    n = HG_params{k}(1);
    m = HG_params{k}(2);
    HG_theory = abs(HG_mode(n, m, X, Y)).^2;

    % Phase-only Hologram
    phase_noise = exp(1i * rand(size(HG_theory))); % random phase noise
    HG_phase_only = abs(ifft2(fft2(HG_theory) .* phase_noise)).^2;

    % CAM Hologram
    amp_mod = 1 + 0.1 * rand(size(HG_theory)); % amplitude modulation
    HG_CAM = abs(ifft2(fft2(HG_theory) .* amp_mod)).^2;

    % save results
    theory_profiles{k + 4} = HG_theory / max(HG_theory(:));
    phase_only_profiles{k + 4} = HG_phase_only / max(HG_phase_only(:));
    cam_profiles{k + 4} = HG_CAM / max(HG_CAM(:));
end

% plot results
figure;
for k = 1:8
    % compute correlation
    C_phase_only = compute_correlation(theory_profiles{k}, phase_only_profiles{k});
    C_CAM = compute_correlation(theory_profiles{k}, cam_profiles{k});

    % Theoretical image
    subplot(3, 8, k);
    imagesc(theory_profiles{k}); axis off; axis equal; colormap parula;
    title(['Theory ', modes{k}]);

    % Phase-only image
    subplot(3, 8, k + 8);
    imagesc(phase_only_profiles{k}); axis off; axis equal; colormap parula;
    title(['Phase-only ', modes{k}, ', C=', num2str(C_phase_only, '%.2f')]);

    % CAM image
    subplot(3, 8, k + 16);
    imagesc(cam_profiles{k}); axis off; axis equal; colormap parula;
    title(['CAM ', modes{k}, ', C=', num2str(C_CAM, '%.2f')]);
end

% custom correlation computation function
function C = compute_correlation(A, B)
    A_mean = mean(A(:));
    B_mean = mean(B(:));
    numerator = sum(sum((A - A_mean) .* (B - B_mean)));
    denominator = sqrt(sum(sum((A - A_mean).^2)) * sum(sum((B - B_mean).^2)));
    C = numerator / denominator;
end
