% define constants
omega0 = 1; % beam waist radius
zR = 10;    % Rayleigh length

% define auxiliary functions
omega = @(z) omega0 * sqrt(1 + (z / zR)^2); % beam width
R = @(z) z * (1 + (zR / z)^2);              % radius of curvature
gouy = @(z) atan(z / zR);                   % Gouy phase

% define LG mode function (full formula)
LG_mode = @(p, l, rho, phi, z) ...
    sqrt(2 * factorial(p) / (pi * factorial(abs(l) + p))) .* ...
    (sqrt(2) * rho ./ omega(z)).^abs(l) .* ...
    laguerreL(p, abs(l), (2 * rho.^2) ./ omega(z).^2) .* ...
    exp(-rho.^2 ./ omega(z).^2) .* ...
    exp(1i * (2 * p + abs(l) + 1) * gouy(z)) .* ...
    exp(-1i * l * phi);

% define HG mode function (full formula)
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

% plot LG modes
figure;
p_values = [0, 2, 3]; % p values
l_values = [0, 1, 2, 4]; % l values
for row = 1:length(p_values)
    for col = 1:length(l_values)
        p = p_values(row);
        l = l_values(col);
        intensity = abs(LG_mode(p, l, rho, phi, z)).^2;
        subplot(length(p_values), length(l_values), (row - 1) * length(l_values) + col);
        imagesc(x, y, intensity);
        axis off;
        axis square;
        % sgtitle（l）
        if row == 1
            title(['$l=', num2str(l), '$'], 'Interpreter', 'latex');
        end
        % note on the left side of each row
        if col == 1
            text(-0.5, 0.5, ['$p=', num2str(p), '$'],  'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Interpreter', 'latex');
        end
    end
end
sgtitle('LG Modes Intensity Profiles');
colormap parula;

% plot HG modes
figure;
m_values = [0, 1, 2]; % m values
n_values = [0, 1, 2, 3]; % n values

for row = 1:length(m_values)
    for col = 1:length(n_values)
        m = m_values(row);
        n = n_values(col);
        intensity = abs(HG_mode(n, m, X, Y, z)).^2;
        subplot(length(m_values), length(n_values), (row - 1) * length(n_values) + col);
        imagesc(x, y, intensity);
        axis off;
        axis square;
        % set top title (n)
        if row == 1
            title(['$n=', num2str(n), '$'], 'Interpreter', 'latex');
        end
        % note on the left side of each row
        if col == 1
            text(-0.5, 0.5, ['$m=', num2str(m), '$'],  'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Interpreter', 'latex');
        end
    end
end
sgtitle('HG Modes Intensity Profiles', 'Interpreter', 'latex');
colormap parula;
