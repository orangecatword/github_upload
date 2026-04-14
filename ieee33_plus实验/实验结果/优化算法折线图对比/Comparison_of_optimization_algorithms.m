% 1. 定义数据
iterations = 1:20; % 迭代次数

pso = [6.40723, 6.40723, 6.40723, 5.39776, 5.04885, 4.87292, 4.59763, 4.35404, 4.28651, 4.28651, 4.28651, 4.28651, 4.28651, 4.28651, 4.28651, 4.28651, 4.25089, 4.25089, 4.25089, 4.25089];
gwo = [4.49170799, 4.49170799, 4.162645554, 3.78662325, 3.763296574, 3.616345945, 3.616345945, 3.616345945, 3.616345945, 3.579057303, 3.579057303, 3.579057303, 3.579057303, 3.579057303, 3.579057303, 3.579057303, 3.579057303, 3.571030545, 3.571030545, 3.571030545];

% 2. 创建图形窗口
figure('Color', [1 1 1]); 
hold on;
grid on;

% 3. 绘制直线图 
% 注意：为了避免标记点(Marker)太挤，每隔2个点绘制一个Marker ('MarkerIndices')
m_step = 1; % 如果点太密可以改为2
plot(iterations, pso, '-o', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));
plot(iterations, gwo, '-s', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));

% 4. 设置坐标轴标签与刻度
xlabel('iterations', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Objective Function', 'FontSize', 12, 'FontWeight', 'bold');
xticks(iterations(1:2:end)); % 隔一个显示一个数字，防止横轴数字太挤

% 【改进点1】显式设置坐标轴字体不为斜体，并设置字体为 Arial 或 Helvetica
set(gca, 'FontAngle', 'normal', 'FontName', 'Arial', 'FontSize', 10);

% 【改进点2】使得间隔明显
% 由于下方三条线在1.5附近非常接近，可以通过设置纵轴范围或手动缩放
% 这里我们设置纵轴范围略大于数据最大值，并开启次要网格线增加层次感
ylim([1.4, 7.0]); 
set(gca, 'YMinorGrid', 'on'); 

% 5. 添加图例
legend({'PSO', 'GWO'}, ...
       'Location', 'northeast', 'FontSize', 10, 'TextColor', 'black');

% 6. 图形美化
set(gca, 'Box', 'on', 'LineWidth', 1.2); 
hold off;