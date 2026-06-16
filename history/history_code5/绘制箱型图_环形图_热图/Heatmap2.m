% 创建数据表
students = {'Alice', 'Bob', 'Charlie', 'Diana'};
subjects = {'Math', 'Physics', 'Chemistry', 'Biology', 'History'};
grades = [90, 85, 78, 92, 88; 
          76, 92, 85, 79, 82;
          88, 79, 91, 84, 76;
          95, 87, 83, 90, 81];

% 创建热力图
figure
h = heatmap(subjects, students, grades);

% 自定义设置
h.Title = 'Student Performance by Subject';
h.Colormap = flipud(hot); % 使用hot颜色映射并反转
h.ColorLimits = [70, 100]; % 固定颜色范围
h.CellLabelFormat = '%d%%'; % 显示百分比格式
h.FontSize = 12;
h.MissingDataColor = [0.8 0.8 0.8]; % 缺失值颜色
h.MissingDataLabel = 'Absent'; % 缺失值标签

% 添加注释
annotation('textbox', [0.1, 0.01, 0.8, 0.05],...
           'String', 'Data Source: Final Exam Records',...
           'EdgeColor', 'none',...
           'HorizontalAlignment', 'center');