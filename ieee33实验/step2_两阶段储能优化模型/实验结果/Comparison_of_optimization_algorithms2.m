% 1. 定义数据
iterations = 1:20; % 迭代次数

bpso_pso = [1.5443624, 1.5443624, 1.5443624, 1.5443624, 1.543158461, 1.542989431, 1.542989431, 1.542989431, 1.542989431, 1.542989431, 1.504701519, 1.504701519, 1.504701519, 1.504701519, 1.504701519, 1.504701519, 1.504701519, 1.504701519, 1.504701519, 1.504701519];
bpso_gwo = [1.5443624, 1.542895584, 1.542895584, 1.533430395, 1.533430395, 1.533430395, 1.533430395, 1.509888383, 1.509888383, 1.509888383, 1.503811036, 1.503811036, 1.503811036, 1.503811036, 1.50183051, 1.50183051, 1.50183051, 1.50183051, 1.50183051, 1.50183051];
bpso_egwo = [1.5443624, 1.5443624, 1.5443624, 1.540208128, 1.533795389, 1.533795389, 1.533795389, 1.512081417, 1.512081417, 1.500872535, 1.500872535, 1.500872535, 1.500872535, 1.500872535, 1.494710256, 1.494710256, 1.494710256, 1.494710256, 1.494710256, 1.494710256];

% 2. 创建图形窗口
figure('Color', [1 1 1]); 
hold on;
grid on;

% 3. 绘制直线图 
% 注意：为了避免标记点(Marker)太挤，每隔2个点绘制一个Marker ('MarkerIndices')
m_step = 1; % 如果点太密可以改为2
plot(iterations, bpso_pso, '-^', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));
plot(iterations, bpso_gwo, '-d', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));
plot(iterations, bpso_egwo, '-v', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));

% 4. 设置坐标轴标签与刻度
xlabel('iterations', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Objective Function', 'FontSize', 12, 'FontWeight', 'bold');
xticks(iterations(1:2:end)); % 隔一个显示一个数字，防止横轴数字太挤

% 【改进点1】显式设置坐标轴字体不为斜体，并设置字体为 Arial 或 Helvetica
set(gca, 'FontAngle', 'normal', 'FontName', 'Arial', 'FontSize', 10);

% 【改进点2】使得间隔明显
% 由于下方三条线在1.5附近非常接近，可以通过设置纵轴范围或手动缩放
% 这里我们设置纵轴范围略大于数据最大值，并开启次要网格线增加层次感
ylim([1.49, 1.55]); 
set(gca, 'YMinorGrid', 'on'); 

% 5. 添加图例
legend({'BPSO+PSO', 'BPSO+GWO','BPSO+EGWO'}, ...
       'Location', 'northeast', 'FontSize', 10, 'TextColor', 'black');

% 6. 图形美化
set(gca, 'Box', 'on', 'LineWidth', 1.2); 
hold off;