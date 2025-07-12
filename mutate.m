%% 变异操作
function offspring = mutate(offspring, pm, var_min, var_max)
    % 输入:
    % offspring: 交叉后的子代种群
    % pm: 变异概率
    % var_min, var_max: 变量的上下界
    % 输出:
    % offspring: 变异后的子代

    n = size(offspring, 1); % 子代种群大小
    var_num = size(offspring, 2); % 变量数量
    eta_m = 20; % 变异的分布指数

    for i = 1:n
        for j = 1:var_num
            if rand <= pm
                % 计算变异扰动
                delta = (rand() * 2 - 1) * (var_max(j) - var_min(j)) / (eta_m + 1);
                offspring(i, j) = offspring(i, j) + delta;

                % 限制在变量上下界内
                offspring(i, j) = max(offspring(i, j), var_min(j));
                offspring(i, j) = min(offspring(i, j), var_max(j));
            end
        end
        % 检查并调整子代以满足限制条件
        offspring(i, :) = enforceConstraints(offspring(i, :));
    end
end

function individual = enforceConstraints(individual)
    % 检查并调整个体以满足限制条件
    % 限制条件:第一到第31列元素与第32到第62列元素间一一对应,要么都为0元素,要么都为非零元素

    % 检查每一列是否满足条件
    for k = 1:31
        if (individual(k) == 0 && individual(k + 31) ~= 0) || (individual(k) ~= 0 && individual(k + 31) == 0)
            % 如果不满足条件,调整其中一个元素使其满足条件
            if individual(k) == 0
                individual(k + 31) = 0;
            else
                individual(k) = 0;
            end
        end
    end
end
