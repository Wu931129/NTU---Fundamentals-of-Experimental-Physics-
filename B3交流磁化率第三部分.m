% --- 讀取資料 ---
data = readmatrix('1008 Xray.txt');
x = data(:,1);   % 角度 (2θ)
y = data(:,2);   % 強度

% --- 限制分析區間 (33° ~ 38°) ---
mask = (x >= 33) & (x <= 38);
x_fit = x(mask);
y_fit = y(mask);

% --- 繪出原始資料 ---
figure;
plot(x, y, 'b'); hold on;
plot(x_fit, y_fit, 'r', 'LineWidth', 1.2);
xlabel('2\theta (degree)');
ylabel('Intensity (a.u.)');
title('XRD Peak (33°–38°)');
grid on;

% --- Gaussian 擬合 ---
% Gaussian model: y = a1*exp(-((x-b1)/c1)^2) + d1
gaussEqn = 'a1*exp(-((x-b1)/c1)^2) + d1';
startPoints = [max(y_fit)-min(y_fit), x_fit(y_fit==max(y_fit)), 0.1, min(y_fit)]; % 初始值
f = fit(x_fit, y_fit, gaussEqn, 'Start', startPoints);

% --- 顯示擬合結果 ---
plot(f, x_fit, y_fit);
legend('原始資料','擬合區間','Gaussian 擬合');

% --- 計算半高寬 (FWHM) ---
% Gaussian form: y = a * exp(-((x-b)/c)^2)
% => FWHM = 2*sqrt(ln(2)) * c
FWHM = 2*sqrt(log(2)) * f.c1;

% --- 顯示結果 ---
fprintf('擬合中心位置 (2θ) = %.4f°\n', f.b1);
fprintf('半高寬 (FWHM) = %.4f°\n', FWHM);
fprintf('峰值強度 = %.4f\n', f.a1 + f.d1);

% --- 顯示擬合曲線疊加圖 ---
figure;
plot(x_fit, y_fit, 'bo'); hold on;
plot(f, 'r-');
xlabel('2\theta (degree)');
ylabel('Intensity (a.u.)');
title(sprintf('Gaussian Fit: Center = %.3f°, FWHM = %.3f°', f.b1, FWHM));
legend('Data','Gaussian fit');
grid on;
