%% 交叉操作
function offspring = crossover(parents, pc, var_num, var_min, var_max)  
    % 输入:
    % parents: 父代种群
    % pc: 交叉概率
    % var_num: 变量数量
    % var_min, var_max: 变量的上下界
    % 输出:
    % offspring: 交叉生成的子代
    parents
    n = size(parents, 1); % 父代种群大小
    offspring = zeros(n, var_num); % 初始化子代种群
    eta_c = 20; % SBX 的分布指数

    for i = 1:2:n-1
        if rand <= pc
            % 随机选择两个父代个体
            parent1 = parents(i, :);
            parent2 = parents(i+1, :);

            % 计算交叉点
            beta = zeros(1, var_num);
            for j = 1:var_num
                if rand <= 0.5
                    beta(j) = (2 * rand())^(1 / (eta_c + 1));
                else
                    beta(j) = (1 / (2 * (1 - rand())))^(1 / (eta_c + 1));
                end
            end

            % 生成子代
            offspring(i, :) = 0.5 * ((1 + beta) .* parent1 + (1 - beta) .* parent2);
            offspring(i+1, :) = 0.5 * ((1 - beta) .* parent1 + (1 + beta) .* parent2);

            % 限制子代在变量上下界内
            offspring(i, :) = max(offspring(i, :), var_min);
            offspring(i, :) = min(offspring(i, :), var_max);
            offspring(i+1, :) = max(offspring(i+1, :), var_min);
            offspring(i+1, :) = min(offspring(i+1, :), var_max);
            % 检查并调整子代以满足限制条件
            offspring(i, :) = enforceConstraints(offspring(i, :));
            offspring(i+1, :) = enforceConstraints(offspring(i+1, :));
        else
            % 如果不进行交叉,直接复制父代
            offspring(i, :) = parents(i, :);
            offspring(i+1, :) = parents(i+1, :);
        end
    end
end
