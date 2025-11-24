%% 拥挤度计算
function crowding_dist = calculate_crowding_distance(obj, fronts)
    % obj: 目标函数值矩阵
    % fronts: Pareto 前沿集合

    n = size(obj, 1); % 种群大小
    m = size(obj, 2); % 目标函数数量  
    crowding_dist = zeros(n, 1); % 初始化拥挤度距离向量

    for i = 1:length(fronts)
        current_front = fronts{i};
        f_min = min(obj(current_front, :), [], 1);
        f_max = max(obj(current_front, :), [], 1);
        for j = 1:m
            % 获取排序后的索引(关键修正)
            [~, sorted_indices] = sortrows(obj(current_front, j));
            sorted_front = current_front(sorted_indices);

            % 处理首尾个体
            crowding_dist(sorted_front(1)) = inf;
            crowding_dist(sorted_front(end)) = inf;

            % 归一化分母处理
            if (f_max(j) - f_min(j)) > 0
                normalization = f_max(j) - f_min(j);
            else
                normalization = 1;
            end

            % 计算中间个体的拥挤度
            for k = 2:(length(sorted_front) - 1)
                crowding_dist(sorted_front(k)) = crowding_dist(sorted_front(k)) + ...
                    (obj(sorted_front(k+1), j) - obj(sorted_front(k-1), j)) / normalization;
            end
        end
    end
end
