PSO = [4.25089	1.78956	2.61619	2.61619	3.78038	1.86393	2.36607	3.71532	4.16111	1.70745	2.06679	5.01526	1.85085	2.47737	2.28065	1.79148	2.36331	1.62517	4.5471	1.99269];
GWO = [3.571030545	1.51232	1.78019	1.67935	1.75107	3.52187	1.67976	4.04163	1.51282	1.77109	1.48337	3.68746	3.6819	1.62721	1.51268	1.51297	3.57286	2.36003	1.60919	1.5298];   
  
% 合并数据
salesData = {PSO, GWO};
groupNames = {'PSO', 'GWO'};
% 创建位置向量
positions = [1, 2]; % 自定义位置避免重叠
% 绘制箱型图
figure
hold on
for i = 1:2
    boxplot(salesData{i}, 'Positions', positions(i), ...
        'Widths', 0.4, 'Colors', rand(1,3))
end

% 美化图形
set(gca, 'XTick', positions, 'XTickLabel', groupNames)
% title('各个算法对目标函数的优化结果')
% ylabel('Objective Function')
ylabel('Objective Function', 'FontWeight', 'bold', 'FontSize', 12)
% xlabel('优化算法')
grid on

% 添加数据点
for i = 1:2
    % 添加抖动点避免重叠
    jitter = 0.05*randn(size(salesData{i}));
    plot(positions(i) + jitter, salesData{i}, 'o', ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], ...
        'MarkerEdgeColor', 'none')
end
hold off