%% 3维PSO算法
function[fy_max] = IEEE123simplePSO_3(~,~)

%% 参数设置
N = 10;                 % 种群个数
d = 3;                 % 维度
ger = 10;              % 最大迭代次数
limit = [0, 10];  % 位置边界
vlimit = [-1, 1];   % 速度边界
c1 = 0.5;               % 自我学习因子
c2 = 0.5;               % 群体学习因子
w = 0.8;                % 惯性权重
%% 新增自适应权重
% 定义权重范围
w_max = 1; 
% w_min = 0.4;

%% 初始化种群
x = limit(1) + (limit(2) - limit(1)) .*rand(N, d);
v = 0.1*(rand(N, d)-0.5);
x_max = x;                     % 个体历史最优位置
y_max = zeros(1, d);           % 全局最优位置
fx_max = ones(N, 1)*inf;       % 个体历史最优适应度
fy_max = inf;                  % 全局最优适应度
record = zeros(ger,1);      % 记录适应度变化

%% PSO迭代
for iter = 1:ger
    %% --- 方案 A：线性递减 (简单常用) ---
    % w = w_max - (w_max - w_min) * (iter / ger);
    tic
    %% 外层目标函数
    fx = down_capacity(x);
    toc
    % 更新个体最优
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better, :) = x(better, :);
    % 更新全局最优
    [minval, nmin] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max  = x_max(nmin, :);
    end 

    % 速度与位置更新
    v = w*v + c1*rand*(x_max - x) + c2*rand*(repmat(y_max, N, 1) - x);
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    x = x + v;
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    % 记录与可视化
    figure(5)
    record(iter) = fy_max;
    subplot(1,2,1);
    bar(y_max);
    title(['迭代', num2str(iter), ' - 全局最优非零元素分布']);
    xlabel('变量索引'); ylabel('值');
    subplot(1,2,2);
    plot(record(1:iter), 'LineWidth', 1.5);
    title('最优适应度进化过程');
    xlabel('迭代次数'); ylabel('最优适应度');
    pause(0.05);

    %% --- 方案 C：基于成功率的自适应 (进阶) ---
    % 如果这一代全局最优 fy_max 没变，说明陷入停滞，适当增大 w 鼓励跳出
    if iter > 1 && record(iter) == record(iter-1)
        w = min(w + 0.05, w_max);
    end
end

%% 输出结果
disp(['最优值：', num2str(fy_max)]);
disp(['最优变量（3维显示）：', num2str(y_max(1:3))]);
disp(['非零变量索引：', num2str(find(y_max~=0))]);
disp(['非零变量取值：', num2str(y_max(y_max~=0))]);

figure(6);
subplot(1,2,1); bar(y_max); title('最终最优变量分布');
subplot(1,2,2); plot(record, 'LineWidth', 1.5); title('PSO适应度收敛曲线');
end
