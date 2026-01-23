%% 量子粒子群算法（定容）
%% 33维QPSO算法（每个个体仅3个非零变量）
function[fy_max] = simpleQPSO(~,~)
%% 参数设置
N = 20;                 % 种群个数
d = 3;                 % 维度
ger = 10;              % 最大迭代次数 (为观察收敛，将迭代次数增大)
limit = [0.8, 1.2];  % 位置边界
% **改进 1：移除速度相关参数和边界**
% vlimit = [-0.02, 0.02];   % [PSO] 速度边界，QPSO移除
% w = 0.8;                % [PSO] 惯性权重，QPSO移除
% c1 = 0.5;               % [PSO] 自我学习因子，QPSO移除
% c2 = 0.5;               % [PSO] 群体学习因子，QPSO移除

% **改进 2：引入 QPSO 关键参数**
c_max = 1.0;            % 收缩-膨胀因子 C 的最大值
c_min = 0.5;            % 收缩-膨胀因子 C 的最小值

%% 初始化种群
x = limit(1) + (limit(2) - limit(1)) .*rand(N, d);
% v = 0.04*rand(N, d); % [PSO] 速度初始化，QPSO移除
x_max = x;                     % 个体历史最优位置
y_max = zeros(1, d);           % 全局最优位置
fx_max = ones(N, 1)*inf;       % 个体历史最优适应度
fy_max = inf;                  % 全局最优适应度
record = zeros(ger,1);      % 记录适应度变化

%% QPSO迭代
tic
for iter = 1:ger
    % 目标函数计算 (保留不变)
    fx = down_capacity(x);
    % 更新个体最优 (保留不变)
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better, :) = x(better, :);
    % 更新全局最优 (保留不变)
    [minval, nmin] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max = x_max(nmin, :);
    end
    
    % **改进 3：计算平均最佳位置 (m_best)**
    m_best = mean(x_max); 
    
    % **改进 4：计算收缩-膨胀因子 C (线性递减策略)**
    % QPSO 关键参数，平衡全局探索和局部开发
    c = c_max - (c_max - c_min) * (iter / ger); 
    
    % **改进 5：QPSO 位置更新 (核心改变)**
    x_new = zeros(N, d);
    
    for i = 1:N
        % 计算吸引点 P_i
        phi = rand(1, d); % [0, 1] 上的随机数
        P_i = phi .* x_max(i, :) + (1 - phi) .* y_max;
        
        % 计算 L_i
        L_i = 2 * c * abs(m_best - x(i, :));
        
        % 计算 u 和 ln(1/u)
        u = rand(1, d); 
        ln_term = log(1./u); 
        
        % QPSO 位置更新公式 (加入随机方向 +/-)
        % 随机决定正负号 (50% 概率)
        sign_term = sign(rand(1, d) - 0.5); 
        
        % 更新位置
        x_new(i, :) = P_i + sign_term .* (L_i / 2) .* ln_term;
    end
    x = x_new; % 用新位置替换旧位置
    
    % 边界处理 (保留不变)
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
 
    % 记录与可视化 (保留不变)
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
toc
%% 输出结果 (保留不变)
disp(['最优值：', num2str(fy_max)]);
disp(['最优变量（仅前10维显示）：', num2str(y_max(1:10))]);
disp(['非零变量索引：', num2str(find(y_max~=0))]);
disp(['非零变量取值：', num2str(y_max(y_max~=0))]);
figure;
subplot(1,2,1); bar(y_max); title('最终最优变量分布');
subplot(1,2,2); plot(record, 'LineWidth', 1.5); title('适应度收敛曲线');