%% 非支配排序
    function [fronts, rank] = non_dominated_sort(obj)

    % obj: 目标函数值矩阵,每一行对应一个个体,每一列对应一个目标函数
    n = size(obj, 1); % 种群大小
    fronts = {}; % 初始化 Pareto 前沿集合
    rank = zeros(n, 1); % 初始化等级向量

    % 初始化每个个体的支配计数和被支配个体集合
    domination_count = zeros(n, 1);
    dominated_sets = cell(n, 1);

    % 计算每个个体的支配计数和被支配个体集合
    for i = 1:n
        for j = 1:n
            if i ~= j
                if all(obj(i, :) <= obj(j, :)) && any(obj(i, :) < obj(j, :))
                    domination_count(j) = domination_count(j) + 1;
                    dominated_sets{i} = [dominated_sets{i}, j];
                end
            end
        end
    end

    % 根据支配计数确定 Pareto 前沿
    front_index = 1;
    remaining = true(n, 1);  % 标记未处理的个体

    while any(remaining)
        % 找到当前前沿的个体(domination_count为0且未被处理)
        current_front = find(domination_count == 0 & remaining);
        if isempty(current_front)
            break;
        end
        fronts{front_index} = current_front;
        
        % 处理当前前沿中的每个个体
        for i = current_front
            rank(i) = front_index;
            remaining(i) = false;  % 标记为已处理
            
            % 更新被支配个体的计数器
            current_set = dominated_sets{i};
            if ~isempty(current_set)
            % if ~isempty(dominated_sets{i})
                current_set = dominated_sets{i};
                for j = 1:length(current_set)
                    domination_count(current_set(j)) = domination_count(current_set(j)) - 1;
                end
            end
        end
        front_index = front_index + 1;
    end
end

