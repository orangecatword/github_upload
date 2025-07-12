%% 选择操作
    function parents = select_parents(population, rank, crowding_dist)
    % 输入:
    % population: 当前种群,每一行是一个个体
    % rank: 非支配等级向量
    % crowding_dist: 拥挤度距离向量
    % 输出:
    % parents: 选择的父代个体

    % 想法:有20套解,但每套解包括pop和Ees_start无直接对应关系,先将其组合起来在进行后续操作

    n = size(population, 1); % 种群大小
    parents = zeros(n, size(population, 2)); % 初始化父代种群

    % 按非支配等级和拥挤度距离排序
    [~, sorted_indices] = sort(rank); % 按非支配等级排序
    sorted_pop = population(sorted_indices, :);
    sorted_rank = rank(sorted_indices);
    sorted_crowding_dist = crowding_dist(sorted_indices);

    % 优先选择非支配等级低的个体
    current_rank = 1;
    index = 1;
    while index <= n
        % 找到当前等级的所有个体
        current_indices = find(sorted_rank == current_rank);
        if ~isempty(current_indices)
            % 如果当前等级的个体数超过剩余需要选择的个体数
            if index + length(current_indices) > n
                % 按拥挤度距离降序选择剩余的个体
                [~, sorted_crowding_indices] = sort(sorted_crowding_dist(current_indices), 'descend');
                selected_indices = current_indices(sorted_crowding_indices(1:(n - index + 1)));
            else
                selected_indices = current_indices;
            end
            % 将选择的个体加入父代
            parents(index:index + length(selected_indices) - 1, :) = sorted_pop(selected_indices, :);
            index = index + length(selected_indices);
        end
        current_rank = current_rank + 1;
    end
end
