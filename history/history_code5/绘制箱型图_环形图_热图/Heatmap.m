% 创建5×5随机矩阵 (0~1)
data = rand(5);

% 绘制基础热力图
figure
h = heatmap(data);

% 添加标题和标签
h.Title = 'Basic Random Matrix Heatmap';
h.XLabel = 'Columns';
h.YLabel = 'Rows';

% 自定义颜色
h.Colormap = parula; % 使用parula颜色映射
colorbar % 显示颜色条