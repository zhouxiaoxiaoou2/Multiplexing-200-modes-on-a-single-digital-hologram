clc; clear; close all;
rng(1); % set random seed for reproducibility

% ========== Define grid ==========
Nx = 128; Ny = 128; % resolution of each mode
x = linspace(-5, 5, Nx);
y = linspace(-5, 5, Ny);
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2); % polar coordinate rho
phi = atan2(Y, X); % polar coordinate phi

% ========== Physical parameters ==========
omega0 = 2; % **beam waist radius**
w0=1;

% ========== Define HG and LG modes ==========
HG_mode = @(m, n, X, Y, w0) ...
    (1 / sqrt(2^m * factorial(m) * 2^n * factorial(n) * pi)) .* ...
    exp(-(X.^2 + Y.^2) / w0^2) .* ...
    hermiteH(n, sqrt(2) * X / w0) .* ...
    hermiteH(m, sqrt(2) * Y / w0);

LG_mode = @(p, l, rho, phi) ...
    (1 / omega0) * sqrt((2 * factorial(p)) / (pi * factorial(abs(l) + p))) .* ...
    (sqrt(2) * rho / omega0).^abs(l) .* ...
    laguerreL(p, abs(l), (2 * rho.^2) / omega0^2) .* ...
    exp(-rho.^2 / omega0^2) .* exp(-1i * l * phi);

% ========== Set 42 mixed modes ==========
num_modes = 42;
grid_size_x = 6; % **6 rows**
grid_size_y = 7; % **7 columns**
spacing = 20; % adjust spacing between modes

% **Randomly decide mode types** (half LG, half HG)
mode_types = randi([0, 1], num_modes, 1); % 0: HG, 1: LG

% **HG mode random parameters**
m_values = randi([0, 2], num_modes, 1);
n_values = randi([0, 2], num_modes, 1);

% **LG mode random parameters**
p_values = randi([0, 2], num_modes, 1);
l_values = randi([-2, 2], num_modes, 1); % **allow negative values for symmetry**

% ========== Generate 42 modes ==========
mixed_transformed = cell(num_modes, 1);
for i = 1:num_modes
    % Choose HG or LG
    if mode_types(i) == 0
        amplitude = HG_mode(m_values(i), n_values(i), X, Y, omega0); % HG mode
    else
        amplitude = LG_mode(p_values(i), l_values(i), rho, phi); % LG mode
    end

    % Store mode
    mixed_transformed{i} = abs(amplitude);
end

% ========== Compute total size needed for arrangement ==========
full_Nx = grid_size_y * (Nx + spacing);
full_Ny = grid_size_x * (Ny + spacing);
final_image = zeros(full_Ny, full_Nx); % **Fix dimensions**

% ========== Fix mode arrangement (6×7 grid) ==========
random_indices = randperm(num_modes); % **Shuffle mode order**
for i = 1:num_modes
    row = mod(i-1, grid_size_x) + 1; % Compute row index
    col = floor((i-1) / grid_size_x) + 1; % Compute column index

    % Compute insertion area
    x_start = (col-1) * (Nx + spacing) + 1;
    y_start = (row-1) * (Ny + spacing) + 1;
    x_range = x_start:(x_start + Nx - 1);
    y_range = y_start:(y_start + Ny - 1);

    % Insert into grid
    final_image(y_range, x_range) = mixed_transformed{random_indices(i)};
end

% ========== Display final result ==========
figure;
imagesc(final_image);
colormap parula; colorbar;
axis off; axis image;
title('42 Mixed LG + HG Modes in 6×7 Grid (Fixed Row/Col Indexing)');
