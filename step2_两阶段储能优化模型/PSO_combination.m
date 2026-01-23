%% 组合编码 PSO
% 不是让每个解在33维里面挑3个为1;而是一个解三个索引,索引在1~33的范围游荡
function [fy_max] = PSO_combination(~,~)

%% ================= 参数设置 =================
N   = 10;          % 粒子数
ger = 10;          % 迭代次数
k   = 3;           % 每个粒子选择 3 个位置
D   = 33;          % 候选位置总数

w  = 0.7;          % 惯性权重
c1 = 1.5;          % 个体学习因子
c2 = 1.5;          % 群体学习因子

vlimit = [-5,5];   % 速度边界

%% ================= 初始化 =====================
x = zeros(N,k);    % 粒子位置（索引编码）
v = zeros(N,k);    % 粒子速度

for i = 1:N
    x(i,:) = sort(randperm(D,k));   % 随机选 3 个不同位置
    v(i,:) = randn(1,k);
end

x_max  = x;
fx_max = inf(N,1);

y_max  = zeros(1,k);
fy_max = inf;

record = zeros(ger,1);

%% ================= PSO 主循环 =================
for iter = 1:ger

    %% --------- 目标函数 ----------
    % 把“索引编码”映射为 0/1 决策向量
    Xbin = zeros(N,D);
    for i = 1:N
        Xbin(i, x(i,:)) = 1;
    end

    fx = down_position(Xbin);

    %% --------- 个体最优 ----------
    better = fx < fx_max;
    fx_max(better) = fx(better);
    x_max(better,:) = x(better,:);

    %% --------- 全局最优 ----------
    [minval, idx] = min(fx_max);
    if minval < fy_max
        fy_max = minval;
        y_max  = x_max(idx,:);
    end

    %% --------- 速度 & 位置更新 ----------
    for i = 1:N
        r1 = rand(1,k);
        r2 = rand(1,k);

        v(i,:) = w*v(i,:) ...
            + c1*r1.*(x_max(i,:) - x(i,:)) ...
            + c2*r2.*(y_max - x(i,:));

        v(i,:) = max(min(v(i,:), vlimit(2)), vlimit(1));

        % 更新位置（连续 → 离散）
        x_new = round(x(i,:) + v(i,:));

        % 边界限制
        x_new(x_new < 1) = 1;
        x_new(x_new > D) = D;

        % 去重 & 修复（组合约束）
        x_new = unique(x_new,'stable');

        while length(x_new) < k
            cand = randi(D);
            if ~ismember(cand,x_new)
                x_new(end+1) = cand;
            end
        end

        x(i,:) = sort(x_new(1:k));
    end

    %% --------- 记录 ----------
    figure(5)
    record(iter) = fy_max;
    subplot(1,2,1);
    bar(accumarray(y_max',1,[D,1]));
    title(['迭代 ',num2str(iter),' - 最优选址分布']);
    xlabel('位置编号'); ylabel('是否被选');

    subplot(1,2,2);
    plot(record(1:iter),'LineWidth',1.5);
    title('适应度收敛');
    xlabel('迭代次数'); ylabel('最优值');
    pause(0.05);
end

%% ================= 输出结果 ==================
disp(['最优目标值：', num2str(fy_max)]);
disp(['最优选址索引：', num2str(y_max)]);

figure(6);
subplot(1,2,1);
bar(accumarray(y_max',1,[D,1]));
title('PSO_combination最终选址结果');

subplot(1,2,2);
plot(record,'LineWidth',1.5);
title('适应度收敛曲线');

end
