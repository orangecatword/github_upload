% 创建示例数据（替换为您的实际数据）
x = linspace(0, 24, 100); % X数据（0-24小时）
y = 10 + 5 * sin(2*pi*x/12) + randn(size(x)); % Y数据（带噪声的周期性数据）

% 将折线图转换为环形图
figure

% 1. 将X数据转换为角度（0到2π范围）
theta = 2 * pi * (x - min(x)) / (max(x) - min(x));

% 2. 将Y数据缩放为半径（反转方向：小值靠近圆心）
r = 1 - (y - min(y)) / (max(y) - min(y)); % 归一化并反转
r = r * 5; % 缩放半径大小（可调整）

% 3. 创建环形图
polarplot(theta, r, 'LineWidth', 2, 'Color', 'b')
hold on

% 4. 添加圆心标记（Y最小值位置）
[min_r, min_idx] = min(r);
polarplot(theta(min_idx), min_r, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')

% 5. 添加坐标标签
ax = gca;
ax.RDir = 'reverse'; % 反转半径方向（数值小靠近圆心）
ax.RAxisLocation = 0; % 半径轴位置（0度处）
ax.ThetaZeroLocation = 'top'; % 0度在顶部（相当于x轴起点）

% 6. 添加标题和标签
title('折线图转环形图', 'FontSize', 14)
text(0, 0, '圆心 (Y最小值)', 'HorizontalAlignment', 'center', 'FontSize', 10)

% 7. 添加参考圆环（辅助解读）
for ref_r = 0:1:max(r)
    polarplot(linspace(0, 2*pi, 100), ref_r*ones(1,100), ':', 'Color', [0.8 0.8 0.8])
end

hold off