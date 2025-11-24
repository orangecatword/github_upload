%% ----------- 简单遗传算法（GA）示例 -----------
clear; clc;

%% ==== 参数设置 ====
popSize = 20;      % 种群规模
maxGen  = 10;      % 最大迭代次数
pc = 0.8;          % 交叉概率
pm = 0.1;          % 变异概率
lb = 0.1; ub = 0.3;   % 搜索范围
nz = 3;            % 只有三个非0元素
dim = 33;          % 维度

%% ==== 初始化种群 ====
pop = zeros(popSize, dim);
for i = 1:popSize
    idx = randperm(dim, nz);     % 随机选3个位置
    pop(i, idx) = lb + rand(1, nz) * (ub - lb);  % 非零元素为随机值
end


%% 适应度函数
% fitnessFunc = @(x) x .* sin(10*pi*x) + 1;

%% ---- 主循环 ----
for gen = 1:maxGen
    % 计算适应度
    fitness = 2/3*down(pop) + 1/3*100*sum(pop,2);
    %% ==== 选择（轮盘赌选择） ====

    invFit = 1 ./ (fitness + 1e-9);      % 防止除零,适应度越小越容易被选中
    prob = invFit / sum(invFit);         % 概率 ∝ 1/fitness
    cumProb = cumsum(prob);                  % 累积概率

    newPop = zeros(popSize,dim);
    for i = 1:popSize
        r = rand;
        idx = find(cumProb >= r, 1, 'first');
        newPop(i,:) = pop(idx,:);                % 被选中的个体
    end
    pop = newPop;

%% ===== 交叉：两个个体仅交换非零元素（保持3个非零） （交换了值没交换位置）=====
for i = 1:2:popSize-1
    if rand < pc
        p1 = pop(i,:);
        p2 = pop(i+1,:);

        % 找出非零元素的位置
        nz1 = find(p1 ~= 0);
        nz2 = find(p2 ~= 0);

        % 随机决定交换一个或两个非零基因
        k = randi([1,2]);

        % 随机选择要交换的位置
        idx1 = nz1(randperm(3, k));   % p1 中的交换位置
        idx2 = nz2(randperm(3, k));   % p2 中的交换位置

        % 交换对应位置的非零值
        temp = p1(idx1);
        p1(idx1) = p2(idx2);
        p2(idx2) = temp;

        % 更新种群
        pop(i,:)   = p1;
        pop(i+1,:) = p2;
    end
end

    %% ==== 变异（随机微扰） ====
    for i = 1:popSize
        if rand < pm
            % 找到当前非零元素位置
            nzIdx = find(pop(i,:) ~= 0);
            % 随机选择一个现有非零置零
            rmIdx = nzIdx(randi(length(nzIdx)));
            pop(i, rmIdx) = 0;
            % 在零位置上随机增加一个新非零
            zeroIdx = find(pop(i,:) == 0);
            addPos = zeroIdx(randi(length(zeroIdx)));
            pop(i, addPos) = lb + rand * (ub - lb);
        end
    end

    % 保持在范围内
    for i = 1:popSize
        for j = 1:dim
            if pop(i,j) ~= 0
                pop(i,j) = max(pop(i,j), lb);
                pop(i,j) = min(pop(i,j), ub);
            end
        end
    end
    %% 输出当前最优（最小适应度）
    [bestFit, bestIdx] = min(fitness);
    fprintf("第 %d 代：最佳 f = %.6f\n", gen, bestFit);
end

%% ===== 最终结果 =====
bestX = pop(bestIdx,:);
bestF = bestFit;

fprintf("\n最终最佳适应度：%.6f\n", bestF);

%% 可视化收敛曲线
figure;
plot(bestHistory,'LineWidth',2);
xlabel('迭代次数');
ylabel('最优适应度值');
title('GA 收敛曲线');
grid on;
% 改进1：是否每次迭代后保留最优值;改进2：交叉部分位置和数值均进行变化