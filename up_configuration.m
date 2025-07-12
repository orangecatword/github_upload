% 上层配置情况(决策变量:DG接入配电网的位置和容量)


function [new_population, new_obj]= up_configuration(C_res,pk,g,b) 

% 参数设置
    pop_size = 20;       % 种群大小(初始化节点数量)
    max_gen = 60;        % 最大迭代次数
    var_num = 62;        % 变量维度
    obj_num = 2;         % 上层目标函数个数.包括年综合成本最小,电压偏差总量最小

    
    % 可调整参数
    pc = 0.9;            % 交叉概率
    pm = 1/var_num;      % 变异概率

% 参数大小设置(初始化变量)
    % LC 代指位置和容量
    LC_wt = sdpvar(pop_size, 33); % 风机对应位置的容量
    LC_pv = sdpvar(pop_size, 33); % 光伏对应位置的容量

    % 光伏和风电出力
    P_wt = sdpvar(pop_size, 25, 33);
    P_pv = sdpvar(pop_size, 25, 33);

    % 分布式电源的有功出力上限(某个场景的某个节点)
    P_DG = sdpvar(pop_size, 25, 33);
    P_DG_max = sdpvar(pop_size, 25, 33);

    % 不同场景下的节点电压
    U = sdpvar(pop_size, 25, 33);
    P = sdpvar(pop_size, 25, 33); % 场景s时流入节点i的有功功率和无功功率
    Q = sdpvar(pop_size, 25, 33);

    % 可转移负载功率
    P_dr = sdpvar(pop_size, 25, 33);

    % 下层目标函数
    f = sdpvar(pop_size, 1);


    

    obj= sdpvar(20, 2); % 不同种群的目标函数值


    
    % 初始化种群
    % 光伏选址定容
    % 风力选址定容
    [LC_wt, LC_pv] = initialize_population(pop_size);

        

    % 下层优化情况
    for i = 1:pop_size
        [P_wt(i,:,:), P_pv(i,:,:), P_DG(i,:,:), P_DG_max(i,:,:), U(i,:,:), P(i,:,:), Q(i,:,:), P_dr(i,:,:), f(i)] = dw_optimum_stand(LC_wt(i,:), LC_pv(i,:), C_res, pk, g, b); % 得到下层输出的参数
    end
 

    % 进化循环
    population = [LC_wt, LC_pv];

    for gen = 1:max_gen
        % 计算目标函数值
        for i = 1:pop_size
            [obj(i)] = evaluate_population(pop_size, LC_wt, LC_pv, P_wt, P_pv, d, P_en, t, P_loss, P_dr);
        end
        % 非支配排序与拥挤度计算
        [fronts, rank] = non_dominated_sort(obj);
        crowding_dist = calculate_crowding_distance(obj, fronts);
        
        % 选择、交叉、变异生成子代
        parents = select_parents(population, rank, crowding_dist);
        
        offspring = crossover(parents, pc, var_num, var_min, var_max);
        offspring = mutate(offspring, pm, var_min, var_max);
        
        % 合并父代与子代
        combined_pop = [parents; offspring];
        pop = offspring(:, 1:31);
        Ees_start = offspring(:, 32:62);

        [new_obj] = evaluate_population(pop_size,Ppv_c,pop,Ees_start,Ppv_big,Qpv_big,Pdec_big,Pes_dc_big,Pes_c_big,f_big,fes_big);
        combined_obj = [obj; new_obj];

        % 环境选择(精英保留)
        [new_population, new_obj] = environmental_selection(combined_pop, combined_obj, pop_size);
        % 更新所有点的位置
        set(h, 'XData', new_obj(:,1), 'YData', new_obj(:,2), 'ZData', new_obj(:,3));
    
        drawnow; % 立即刷新图形
        pause(0.1); % 控制刷新速度
    end
 
% end