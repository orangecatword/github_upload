PSO = [1.780195064	1.952641115	1.817299305	1.911225231	1.38432606	1.76276315	1.435097733	1.517058476];
GWO = [1.783827697	1.783827697	1.435260917	1.609928583	1.79875809	1.497999017];
EGWO = [1.345543907	1.041749368	1.344791551	1.354903563	1.678273293];
% 合并数据
salesData = {PSO, GWO, EGWO};
groupNames = {'PSO','GWO', 'EGWO'};
% 创建位置向量
positions = [1, 2, 3]; % 自定义位置避免重叠
% 绘制箱型图
figure
hold on
for i = 1:3
    boxplot(salesData{i}, 'Positions', positions(i), ...
        'Widths', 0.4, 'Colors', rand(1,3))
end

% —————核心修改部分—————
% 设置Y轴范围：从 1.4 开始，上限根据数据最大值（约1.7）稍作预留
ylim([0.8, 2]); 
% —————————————————————

% 美化图形
set(gca, 'XTick', positions, 'XTickLabel', groupNames)
% title('各个算法对目标函数的优化结果')
% ylabel('Objective Function')
ylabel('Objective Function', 'FontWeight', 'bold', 'FontSize', 12)
% xlabel('优化算法')
grid on

% 添加数据点
for i = 1:3
    % 添加抖动点避免重叠
    jitter = 0.05*randn(size(salesData{i}));
    plot(positions(i) + jitter, salesData{i}, 'o', ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], ...
        'MarkerEdgeColor', 'none')
end
hold off