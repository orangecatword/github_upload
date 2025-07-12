%% 环境选择
function [new_population, new_obj] = environmental_selection(combined_pop, combined_obj, pop_size)
    % 输入:
    % combined_pop: 合并的种群(父代和子代)
    % combined_obj: 合并种群的目标函数值
    % pop_size: 新种群的大小
    % 输出:
    % new_population: 经过环境选择后的新种群
    % new_obj: 新种群的目标函数值

    % 非支配排序
    [fronts, rank] = non_dominated_sort(combined_obj);

    % 拥挤度计算
    crowding_dist = calculate_crowding_distance(combined_obj, fronts);

    % 初始化新种群
    new_population = zeros(pop_size, size(combined_pop, 2));
    new_obj = zeros(pop_size, size(combined_obj, 2));

    % 当前选择的个体数
    count = 0;

    % 按非支配等级选择个体
    for i = 1:length(fronts)
        front = fronts{i};
        if count + length(front) <= pop_size
            % 如果当前前沿的所有个体都能被选入新种群
            new_population(count + 1:count + length(front), :) = combined_pop(front, :);
            new_obj(count + 1:count + length(front), :) = combined_obj(front, :);
            count = count + length(front);
        else
            % 如果当前前沿的个体数超过剩余需要选择的个体数
            % 按拥挤度距离降序选择剩余的个体
            [~, sorted_indices] = sort(crowding_dist(front), 'descend');
            selected_indices = front(sorted_indices(1:(pop_size - count)));
            new_population(count + 1:pop_size, :) = combined_pop(selected_indices, :);
            new_obj(count + 1:pop_size, :) = combined_obj(selected_indices, :);
            count = pop_size;
            break;
        end
    end
end

% end