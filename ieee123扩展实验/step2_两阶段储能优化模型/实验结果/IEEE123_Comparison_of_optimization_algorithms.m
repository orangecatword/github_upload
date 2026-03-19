% 1. 定义数据
iterations = 1:20; % 迭代次数

pso = [21.59649474	21.59649474	21.59649474	21.59649474	19.82159732	17.52866039	16.60563358	16.00115394	15.75660598	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725	14.32254725];
gwo = [16.27951653	16.27951653	15.33212186	15.33212186	15.33212186	15.22136984	14.40220409	14.2964013	14.29513308	14.29513308	14.19705493	14.19705493	14.19705493	14.19705493	14.19571541	14.19571541	14.19571541	14.19571541	14.19571541	14.19550667];
bpso_pso = [14.3425	14.0107	12.913	12.8589	12.8589	12.8589	12.2808	12.0703	11.8682	11.7935	11.7935	11.7106	11.7106	11.7106	11.6793	11.6793	11.6793	11.6793	11.6793	11.6793];
bpso_gwo = [14.34245781	13.60142086	12.12485024	11.66955634	11.66955634	11.66955634	11.66955634	11.66955634	11.66955634	11.65371564	11.65371564	11.65371564	11.64635682	11.64574968	11.64574968	11.64574968	11.64574968	11.64524307	11.64524307	11.64524307];
bpso_egwo = [14.3425	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478	11.6478];

% 2. 创建图形窗口
figure('Color', [1 1 1]); 
hold on;
grid on;

% 3. 绘制直线图 
% 注意：为了避免标记点(Marker)太挤，每隔2个点绘制一个Marker ('MarkerIndices')
m_step = 1; % 如果点太密可以改为2
plot(iterations, pso, '-o', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));
plot(iterations, gwo, '-s', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerIndices', 1:m_step:length(iterations));
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
ylim([11, 22]); 
set(gca, 'YMinorGrid', 'on'); 

% 5. 添加图例
legend({'PSO', 'GWO','BPSO+PSO', 'BPSO+GWO','BPSO+EGWO'}, ...
       'Location', 'northeast', 'FontSize', 10, 'TextColor', 'black');

% 6. 图形美化
set(gca, 'Box', 'on', 'LineWidth', 1.2); 
hold off;