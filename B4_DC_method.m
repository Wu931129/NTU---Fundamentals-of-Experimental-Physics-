% =============================================================
% 程式碼 3: 結合 LSM 與 ODR 的完整分析腳本
%
% 此腳本將分別執行：
% 1. 傳統最小平方法 (LSM) - 忽略 X 軸誤差
% 2. 正交距離回歸 (ODR) - 考量 X 軸與 Y 軸誤差
%
% 將產生兩張獨立的圖形。
% **不需要額外工具箱**
% =============================================================

clear; clc; close all;

% --- 實驗數據 (定義為行向量以便於 ODR 矩陣運算) ---
I_data = [0.200, 0.507, 1.043, 1.503, 1.996, 2.503, 3.070]';
I_err = [0.070, 0.019, 0.070, 0.040, 0.035, 0.036, 0.019]';

BA_data = [0.0010, 0.0007, 0.0168, 0.0321, 0.0353, 0.0418, 0.0483]';
BA_err = [0.0004, 0.0005, 0.0046, 0.0036, 0.0032, 0.0019, 0.0010]';

n = length(I_data); % 數據點數量

% ===================================
% 第一部分：最小平方法 (LSM) 分析
% ===================================

% --- 1a. 使用 polyfit 執行標準 LSM ---
p_lsm = polyfit(I_data, BA_data, 1);
a_lsm = p_lsm(1);
b_lsm = p_lsm(2);

% --- 1b. 計算 LSM 的標準誤差 ---
y_fit_lsm = a_lsm * I_data + b_lsm;
res_lsm = BA_data - y_fit_lsm;
sigma_y_lsm = sqrt(sum(res_lsm.^2) / (n - 2)); % 殘差標準差
Sxx = sum((I_data - mean(I_data)).^2);
sigma_a_lsm = sigma_y_lsm / sqrt(Sxx);
sigma_b_lsm = sigma_y_lsm * sqrt(sum(I_data.^2) / (n * Sxx));

% ===================================
% 第二部分：正交距離回歸 (ODR) 分析
% ===================================

% --- 2a. ODR 迭代計算 ---
X = [I_data, ones(n, 1)]; % 設計矩陣

% 使用 LSM 結果作為初始猜測值
a_odr = a_lsm;
b_odr = b_lsm;

num_iterations = 10; % 迭代 10 次
for i = 1:num_iterations
    % 1. 計算有效權重 (Effective Weights)
    % w_i = 1 / (sigma_y_i^2 + a_i^2 * sigma_x_i^2)
    w = 1.0 ./ (BA_err.^2 + a_odr^2 * I_err.^2);
    
    % 2. 建立權重矩陣
    W = diag(w);
    
    % 3. 執行加權最小平方法 (WLS)
    % 求解 (X' * W * X) * beta = (X' * W * Y)
    beta_odr = (X' * W * X) \ (X' * W * BA_data);
    
    % 4. 更新斜率 a 以進行下一次迭代
    a_odr = beta_odr(1);
    b_odr = beta_odr(2);
end

% --- 2b. 計算 ODR 的標準誤差 ---
res_odr = BA_data - (a_odr * I_data + b_odr);
cov_matrix = inv(X' * W * X); % 使用最後一次迭代的權重 W

% 殘差加權平方和 (Chi-squared)
chi2 = sum(w .* res_odr.^2);
% 約化卡方 (Reduced Chi-squared)
chi2_per_dof = chi2 / (n - 2); 

% 縮放協方差矩陣 (如果 chi2_per_dof > 1，則放大誤差)
cov_matrix_scaled = cov_matrix * chi2_per_dof;

std_errors = sqrt(diag(cov_matrix_scaled));
sigma_a_odr = std_errors(1);
sigma_b_odr = std_errors(2);

% ===================================
% 第三部分：在命令視窗中顯示結果
% ===================================

fprintf('--- 擬合結果比較 ---\n\n');
fprintf('  方法\t\t 斜率 (a)\t\t 截距 (b)\n');
fprintf('---------------------------------------------------\n');
fprintf('  LSM (忽略X誤差)\t %.5f (±%.5f)\n', a_lsm, sigma_a_lsm);
fprintf('  ODR (考量X,Y誤差)\t %.5f (±%.5f)\n\n', a_odr, sigma_a_odr);
fprintf('  方法\t\t 截距 (b)\n');
fprintf('---------------------------------------------------\n');
fprintf('  LSM (忽略X誤差)\t %.5f (±%.5f)\n', b_lsm, sigma_b_lsm);
fprintf('  ODR (考量X,Y誤差)\t %.5f (±%.5f)\n', b_odr, sigma_b_odr);

% ===================================
% 第四部分：繪製圖形
% ===================================

% --- 4a. 繪製 LSM 圖形 (Figure 1) ---
figure(1);
set(gcf, 'Color', [0.1 0.1 0.1]); % 深灰背景
hold on; grid on; box on;

% 誤差棒與數據點
errorbar(I_data, BA_data, BA_err, BA_err, I_err, I_err, 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'w', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'CapSize', 10);

% LSM 回歸線
x_fit_range = linspace(min(I_data)*0.9, max(I_data)*1.05, 100);
plot(x_fit_range, a_lsm * x_fit_range + b_lsm, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

% 標籤與主題
xlabel('電流 (A)', 'FontSize', 14, 'Color', 'w');
ylabel('角度差 B - A (rad)', 'FontSize', 14, 'Color', 'w');
title('電流 vs. 角度(B - A) 最小平方法 (LSM) 分析', 'FontSize', 15, 'Color', 'w');

% 軸設定
set(gca, 'Color', [0.05 0.05 0.05], 'XColor', 'w', 'YColor', 'w', ...
    'GridColor', [0.5 0.5 0.5], 'LineWidth', 1.2, 'FontSize', 12);
ylim([min(BA_data)-0.005, max(BA_data)+0.005]);

% 顯示 LSM 擬合結果文字
text(0.4, max(BA_data)-0.006, ...
    sprintf('LSM 擬合: y = %.4f(±%.4f)x + %.4f(±%.4f)', a_lsm, sigma_a_lsm, b_lsm, sigma_b_lsm), ...
    'Color', 'w', 'FontSize', 12, 'BackgroundColor', [0.15 0.15 0.15], ...
    'EdgeColor', [0.4 0.4 0.4], 'Margin', 5);

legend({'平均值 (含誤差)', 'LSM 線性擬合'}, 'TextColor', 'w', 'Location', 'northwest');
hold off;

% --- 4b. 繪製 ODR 圖形 (Figure 2) ---
figure(2);
set(gcf, 'Color', [0.1 0.1 0.1]); % 深灰背景
hold on; grid on; box on;

% 誤差棒與數據點
errorbar(I_data, BA_data, BA_err, BA_err, I_err, I_err, 'o', ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'w', ...
    'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'CapSize', 10);

% ODR 回歸線
plot(x_fit_range, a_odr * x_fit_range + b_odr, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

% 標籤與主題
xlabel('電流 (A)', 'FontSize', 14, 'Color', 'w');
ylabel('角度差 B - A (rad)', 'FontSize', 14, 'Color', 'w');
title('電流 vs. 角度(B - A) 正交距離回歸 (ODR) 分析', 'FontSize', 15, 'Color', 'w');

% 軸設定
set(gca, 'Color', [0.05 0.05 0.05], 'XColor', 'w', 'YColor', 'w', ...
    'GridColor', [0.5 0.5 0.5], 'LineWidth', 1.2, 'FontSize', 12);
ylim([min(BA_data)-0.005, max(BA_data)+0.005]);

% 顯示 ODR 擬合結果文字
text(0.4, max(BA_data)-0.006, ...
    sprintf('ODR 擬合: y = %.4f(±%.4f)x + %.4f(±%.4f)', a_odr, sigma_a_odr, b_odr, sigma_b_odr), ...
    'Color', 'w', 'FontSize', 12, 'BackgroundColor', [0.15 0.15 0.15], ...
    'EdgeColor', [0.4 0.4 0.4], 'Margin', 5);

legend({'平均值 (含誤差)', 'ODR 線性擬合'}, 'TextColor', 'w', 'Location', 'northwest');
hold off;