% 生成模拟销售数据
east = randn(100,1)*5000 + 25000;    % 东部地区
west = randn(80,1)*8000 + 30000;     % 西部地区
south = randn(120,1)*4000 + 22000;   % 南部地区
north = randn(90,1)*6000 + 28000;    % 北部地区

% 合并数据
salesData = {east, west, south, north};
groupNames = {'东部', '西部', '南部', '北部'};

% 创建位置向量
positions = [1, 1.5, 2.5, 3]; % 自定义位置避免重叠

% 绘制箱型图
figure
hold on
for i = 1:4
    boxplot(salesData{i}, 'Positions', positions(i), ...
        'Widths', 0.4, 'Colors', rand(1,3))
end
hold off

% 美化图形
set(gca, 'XTick', positions, 'XTickLabel', groupNames)
title('各地区季度销售额分布')
ylabel('销售额（元）')
xlabel('地区')
grid on

% 添加数据点
hold on
for i = 1:4
    % 添加抖动点避免重叠
    jitter = 0.05*randn(size(salesData{i}));
    plot(positions(i) + jitter, salesData{i}, 'o', ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], ...
        'MarkerEdgeColor', 'none')
end
hold off