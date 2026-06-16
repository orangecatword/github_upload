
%% 定义参数变量
iter = 20; % 运行优化算法的次数
idx = zeros(iter, 3); % 储能所在位置
vals= zeros(iter, 3); % 储能容量的大小

f_max = zeros(1, iter); % 目标函数最优值
record = zeros(20, iter); % 目标函数曲线
x = zeros(iter, 33); % 最终变量
for i =1:3
    [f_max(:,i),x(i, :),record(:,i)] = simplePSO;
    % [f_max,x,record] = simpleGWO;
    % [f_max,x,record] = simpleBPSO;
    % [f_max,x,record] = simplePSO_3;
    % [f_max,x,record] = simpleGWO_3;
    % [f_max,x,record] = simpleEGWO_3;
    idx(i, :) = find(x(i, :) ~= 0);
    vals(i, :) = x(i,x(i, :) ~= 0);
end
f_max;
