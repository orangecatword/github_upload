function[objective,es_nodes,Ees_max] = function_up_last(~,~,~)

%% 1.设参
% 筛选的5个节点重要性最高的点
candidate_nodes = [24 23 6 8 30];
N_candidate = length(candidate_nodes);

% 最多安装5套，数量由CPLEX决定
N_es_max = 5;
N_es_min = 1;     

% 储能容量上下限，需要根据实际参数修改
Ees_max_lb = 0.05*ones(N_candidate,1);
Ees_max_ub = 1.00*ones(N_candidate,1);

%% 2.设变量
% 每次运行前清理掉上次运行残留的变量
for i = 1:1
  yalmip('clear')

% 储能安装状态
x_es = binvar(N_candidate,1);

% 每个候选节点的额定容量
Ees_max = sdpvar(N_candidate,1);

%% 3.设约束
Constraints = [];    

% 储能数量约束
Constraints = [Constraints,...
    N_es_min <= sum(x_es) <= N_es_max];

% 未安装节点容量必须为0
Constraints = [Constraints,...
    Ees_max >= Ees_max_lb.*x_es,...
    Ees_max <= Ees_max_ub.*x_es];


ob1 = function_down(x_es,Ees_max,11)+function_down(x_es,Ees_max,23)+function_down(x_es,Ees_max,6)+function_down(x_es,Ees_max,3)+function_down(x_es,Ees_max,20);
% ob2 = function_down(x_es,Ees_max,[29-2])+function_down(x_es,Ees_max,[24-14])+function_down(x_es,Ees_max,[22-7])+function_down(x_es,Ees_max,[12-17])+function_down(x_es,Ees_max,[18-26]);
% ob3 = function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11]);
% ob4 = function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11]);
% ob5 = function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11])+function_down(x_es,Ees_max,[11]);

Ees_cost = 105; % 电池配置价格 单位:万美元/10MWh 
w_ess = 0.003; % 储能成本权重-每年可能发生的 1-2 次 电网断电事故;储能预期使用年限：10 – 15 年
X = w_ess*Ees_cost*sum(Ees_max);

objective = ob1+X;

%% 5.设求解器
ops = sdpsettings('solver', 'cplex', 'verbose', 0);
ops.cplex.preprocessing.presolve = 1; % 启用预处理
ops.cplex.workmem = 4096; 
ops.cplex.mip.strategy.search = 1;  % 使用动态搜索
ops.cplex.mip.tolerances.mipgap = 0.0001;  % 放宽最优间隙
ops.cplex.mip.tolerances.absmipgap = 1e-3; % 绝对间隙
ops.cplex.mip.strategy.heuristicfreq = 50; % 调整启发式频率
ops.cplex.threads = 1;        % ⭐ 关键 限制每个 CPLEX 求解器仅用 1 线程
ops.cplex.parallel = 0;       % ⭐ 关键 禁用 CPLEX 内部并行
ops.cplex.nodefileind = 2;    % 启用节点压缩磁盘文件
sol=optimize(Constraints,objective,ops);
end

objective = value(objective);
x_es_value = round(value(x_es));
es_nodes = candidate_nodes(x_es_value > 0.5);
Ees_max_all = value(Ees_max);
Ees_max = Ees_max_all(x_es_value > 0.5);
