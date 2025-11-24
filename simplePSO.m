%% 33维PSO算法（每个个体仅3个非零变量）
function[fy_max] = up(~,~)


%% 参数设置
N = 50;                 % 种群个数
d = 33;                 % 维度
ger = 2;              % 最大迭代次数
limit = [0.1, 0.3];  % 位置边界
vlimit = [-0.02, 0.02];   % 速度边界
w = 0.8;                % 惯性权重
c1 = 0.5;               % 自我学习因子
c2 = 0.5;               % 群体学习因子

%% 初始化种群
x = limit(1) + (limit(2) - limit(1)) .*rand(N, d);
v = 0.04*rand(N, d);

% 位置范围(从-5~5 到0.1~0.3 缩小50倍)和速度大小需要调节

% 强制每个个体仅3个非零变量
for i = 1:N
    idx = randperm(d, 3);
    mask = false(1, d);
    mask(idx) = true;
    x(i, ~mask) = 0;
end

x_max = x;                     % 个体历史最优位置
y_max = zeros(1, d);           % 全局最优位置
fx_max = ones(N, 1)*inf;       % 个体历史最优适应度
fy_max = inf;                  % 全局最优适应度
record = zeros(ger,1);      % 记录适应度变化

%% PSO迭代
for iter = 1:ger
    tic
    % 此时fx应该是外层目标函数-由于论文中没有提及,则定义w_F,ESS;w_F,load;w_F,net分别为1/3
    %% 外层目标函数
    fx = 2/3*down(x) + 1/3*100*sum(x,2); % 加上储能成本,成本为100美元/千瓦时=100万美元/10MWh
    toc
    % 更新个体最优
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better, :) = x(better, :);
    
    % 更新全局最优
    [minval, nmin] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max = x_max(nmin, :);
    end
    
    % 速度与位置更新
    v = w*v + c1*rand*(x_max - x) + c2*rand*(repmat(y_max, N, 1) - x);
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    x = x + v;
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    
    % 保证每个个体只有3个非零元素-矩阵稀疏化,只保证矩阵中最大的三个元素非零-这不就导致每一行的储能位置固定住了;缺陷2:迭代后期是否应该逐渐收敛 
    for i = 1:N
        [~, idx] = sort(abs(x(i,:)), 'descend');
        mask = false(1, d);
        mask(idx(1:3)) = true;
        x(i, ~mask) = 0;
    end
    
    % 记录与可视化
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
end


%% 输出结果
disp(['最优值：', num2str(fy_max)]);
disp(['最优变量（仅前10维显示）：', num2str(y_max(1:10))]);
disp(['非零变量索引：', num2str(find(y_max~=0))]);
disp(['非零变量取值：', num2str(y_max(y_max~=0))]);

figure;
subplot(1,2,1); bar(y_max); title('最终最优变量分布');
subplot(1,2,2); plot(record, 'LineWidth', 1.5); title('适应度收敛曲线');
