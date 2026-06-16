
% 准备数据
x = [0 22.0899	5.8125	2.30349	0.945554	0.658464	0.65019	0.370279];
iterations = 2:8;

figure; % 新建图形窗口

% 1. 绘制蓝色实线和圆点 ('-bo' 表示实线-蓝色-圆点)
% MarkerFaceColor 指定圆点填充颜色
plot(iterations, x(2:8), '-bo', 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'b', ...
    'MarkerSize', 6); 

hold on; % 保持图形，准备叠加标注点

% 2. 在 x=3 的节点标红色星号 (对应 iterations 的第 2 个元素)
% 'rp' 表示红色五角星 (pentagram)，'MarkerSize' 控制星星大小
plot(4, x(4), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

% 3. 设置坐标轴和网格
grid on; % 开启网格线，与图中效果一致
title('SSE with different number of clusters');
xlabel('Number of clusters');
ylabel('SSE');

% 4. 添加图例 (可选)
legend('SSE values', 'Selected Point (k=4)', 'Location', 'northeast');

% 优化坐标轴范围（可选，根据需要调整）
set(gca, 'XTick', 2:8); % 确保横坐标显示整数刻度