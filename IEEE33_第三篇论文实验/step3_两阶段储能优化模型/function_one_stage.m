% clear 
% clc
% warning off
%% 按照文章统一参数
%% 将时间T从1:4变为1
function[objective,r_load,V_bias,R_load] = function_one_stage(~,~,~)

% 筛选的5个节点重要性最高的点
candidate_nodes = [24 23 6 8 30];
N_candidate = length(candidate_nodes);
Result = [23];

% 最多安装5套，数量由CPLEX决定
N_es_max = 5;
N_es_min = 1;     

% 储能容量上下限，需要根据实际参数修改
Ees_max_lb = 0.05*ones(N_candidate,1);
Ees_max_ub = 1.00*ones(N_candidate,1);


%% 1.设参
mpc = case33bw2;
pload = mpc.pload;% 负荷数据
pload_prim = mpc.pload_prim/(1000*10);% 10为基准值 最后得到的负荷为标幺值
qload_prim = mpc.qload_prim/(1000*10);
a = 3.715; % 单时段所有节点有功容量,MW
b = 2.3; % 单时段所有节点无功容量,MW
pload = pload/a;% 得到各个时段与单时段容量的比例系数
qload = pload/b;% 假设有功负荷曲线与无功负荷曲线相同
pload = pload_prim*pload;% 得到33*24的负荷值,每一个时间段每个节点的负荷
qload = qload_prim*qload;
% 选取9-12时作为故障隔离时间
pload = pload(:,9:12);
qload = qload(:,9:12);
branch = mpc.branch;         
r=branch(:,3);         
x=branch(:,4);            
T = 4; %时段数为24小时-改为4小时
nb = 33;%节点数
nl = 37;%支路数
% nc = 5; %联络开关数
upstream=zeros(33,37);%代表流入节点支路
dnstream=zeros(33,37);%代表流出节点支路

for i=1:37
    upstream(branch(i,2),i)=1;
end

for i=1:37
    dnstream(branch(i,1),i)=1;
end

%% 1. 重新定义物理出力限值(Gemini新加入)
% 初始化所有节点出力上限为 0
Pgmax = zeros(nb, T); 
Qgmax = zeros(nb, T);

% (1) 主网节点 (节点 1)
Pgmax(1, :) = 0.09; % 对应 10MW
Qgmax(1, :) = 0.09;

% (2) 定义分布式电源 (DER) 的物理容量
% 建议容量设得略大于风光出力的峰值，以防逆变器限功率
% 根据Case33bw,参照之前的分析：风机峰值约 1.5MW (0.15pu)，光伏峰值约 1MW (0.1pu)
pv_nodes = [7, 27];
wt_nodes = 12;

Pgmax(pv_nodes, :) = 0.03; 
Qgmax(pv_nodes, :) = 0.03; % 给予一定的无功补偿空间
Pgmax(wt_nodes, :) = 0.03; 
Qgmax(wt_nodes, :) = 0.03;
% Pgmax(es_nodes, :) = 0.04; 
% Qgmax(es_nodes, :) = 0.04;

Pgmax(candidate_nodes,:) = 0.04;
Qgmax(candidate_nodes, :) = 0.04;

%% 2.设变量
for i = 1:1
  yalmip('clear')

% 储能安装状态
x_es = binvar(N_candidate,1);

% 每个候选节点的额定容量
Ees_max = sdpvar(N_candidate,1);

U = sdpvar(nb,T);%电压的平方
Iij = sdpvar(nl,T);%电流的平方
Pij = sdpvar(nl,T);%线路有功
Qij = sdpvar(nl,T);%线路无功
Pg = sdpvar(nb,T);%发电机有功
Qg = sdpvar(nb,T);%发电机无功

lamda=sdpvar(33,T,'full');% 负荷损失的状态变量
Zij=binvar(nl,T);%网架结构               

% 储能状态
% Ies_c = binvar(3, T,'full');% 充电状态
% Ies_dc = binvar(3, T,'full');% 放电状态
Ies_c  = binvar(N_candidate,T,'full');
Ies_dc = binvar(N_candidate,T,'full');


% 储能容量-根据论文修改
% Ees_max = Ees_max';

% 储能设备剩余能量
% Ees = sdpvar(3, T+1,'full'); 
Ees = sdpvar(N_candidate,T+1,'full');

% 储能功率
% Pes_c = sdpvar(3, T,'full');% 充电功率
% Pes_dc = sdpvar(3, T,'full');% 放电功率
Pes_c  = sdpvar(N_candidate,T,'full');
Pes_dc = sdpvar(N_candidate,T,'full');

%% 3.设约束
Constraints = [];    

% 储能数量约束
Constraints = [Constraints,...
    N_es_min <= sum(x_es) <= N_es_max];

% 未安装节点容量必须为0
Constraints = [Constraints,...
    Ees_max >= Ees_max_lb.*x_es,...
    Ees_max <= Ees_max_ub.*x_es];


%% 储能约束-可以取消充放电的功率和最好一致的约束，以增强储能供电能力
x_es_T = repmat(x_es,1,T);

Constraints = [Constraints,...
    Ies_c == 0,...
    Ies_dc == x_es_T,...
    Pes_c == 0,...
    Pes_dc >= 0];

% 放电功率不超过0.1倍额定能量容量
Constraints = [Constraints,...
    Pes_dc <= 0.1*repmat(Ees_max,1,T)];

% 初始SOC
Constraints = [Constraints,...
    Ees(:,1) == 0.5*Ees_max];

% SOC上下限
Constraints = [Constraints,...
    Ees >= 0.2*repmat(Ees_max,1,T+1),...
    Ees <= 0.8*repmat(Ees_max,1,T+1)];

% 能量变化方程
Constraints = [Constraints,...
    Ees(:,2:T+1) == Ees(:,1:T) ...
    + 0.9*Pes_c ...
    - 1.1*Pes_dc];

% Constraints = [Constraints;
%     Ies_c == 0;
%     Ies_dc == 1; % 储能状态约束
%     Pes_c == 0;
%     Pes_dc >= 0;
%     Pes_dc <= 0.1*Ies_dc.*(Ees_max*ones(1,T));
%     Ees(:,1) == 0.5*Ees_max;
%     Ees >= 0.2 * Ees_max*ones(1,T+1);
%     Ees <= 0.8 * Ees_max*ones(1,T+1)
%     ];
% Constraints = [Constraints;
%     Ees(:, 2:T+1) == Ees(:, 1:T) + 0.9*Pes_c - 1.1*Pes_dc];

%% 潮流约束
% 节点功率约束
Constraints = [Constraints, pload.*lamda - Pg == upstream*Pij - upstream*(Iij.*(r*ones(1,T))) - dnstream*Pij];
Constraints = [Constraints, qload.*lamda - Qg == upstream*Qij - upstream*(Iij.*(x*ones(1,T))) - dnstream*Qij];%节点注入无功

m = 1.06*1.06 - 0.94*0.94;
M = (ones(nl,T) - Zij)*m;             
Constraints = [Constraints, U(branch(:,1),:) - U(branch(:,2),:) <=  M + 2*(r*ones(1,T)).*Pij + 2*(x*ones(1,T)).*Qij - ((r.^2 + x.^2))*ones(1,T).*Iij];
Constraints = [Constraints, U(branch(:,1),:) - U(branch(:,2),:) >= -M + 2*(r*ones(1,T)).*Pij + 2*(x*ones(1,T)).*Qij - ((r.^2 + x.^2))*ones(1,T).*Iij];

%节点电压约束               
Constraints = [Constraints, U(1, :) == 1.0]; % 首节点电压为1
Constraints = [Constraints, 0.94^2<= U,U <= 1.06^2]; % 电压波动范围约束

%% 商品流约束-问题：出现了环流情况
% 论文观点为,除开故障状况,否则32条主干支路无法随便关断
% pg_st = [1 7 12 27 es_nodes];
% Fij=sdpvar(37,T,'full'); 
% Wj=sdpvar(7,T,'full'); 
% M_2=100;
for t=1:T
    % for k=1:33
    %     if ~ismember(k,pg_st)
    %         node_out=find(branch(:,1)==k);
    %         node_in=find(branch(:,2)==k);
    %         Constraints=[Constraints,sum(Fij(node_out,t))-sum(Fij(node_in,t))==-1];
    %     else
    %         node_out=find(branch(:,1)==k);
    %         node_in=find(branch(:,2)==k);
    %         Constraints=[Constraints,sum(Fij(node_out,t))-sum(Fij(node_in,t))==Wj(find(pg_st==k),t)]; % 更接近于原本的参考程序
    %     end
    % end
    Constraints=[Constraints,Zij(Result,t) == 0];
end
% Constraints=[Constraints,-M_2.*Zij<=Fij<=M_2.*Zij];
% Constraints=[Constraints,-M_2.*(2-Zij)<=Fij<=M_2.*(2-Zij)];
% Constraints=[Constraints,Wj>=1];

Constraints=[Constraints,sum(Zij,1) <= 32];
Constraints=[Constraints,sum(Zij(1:32,:),1) >= 32-size(Result,2)];


%二阶锥约束  
for i = 1:37
    for t = 1:T
        Constraints = [Constraints;
            cone([2*Pij(i,t); 2*Qij(i,t);Iij(i,t) - U(branch(i,1),t)], Iij(i,t) + U(branch(i,1),t))]; 
    end
end

%% 通用约束
%发电机功率约束           
Constraints = [Constraints, 0 <= Pg,Pg <= Pgmax,-Qgmax <= Qg,Qg <= Qgmax];
% Constraints = [Constraints, Pg(es_nodes, :) == Pes_dc];
Constraints = [Constraints,...
    Pg(candidate_nodes,:) == Pes_dc];

% 新增的针对储能未知情况下的无功功率约束
Constraints = [Constraints,...
Qg(candidate_nodes,:) <= 0.04*repmat(x_es,1,T);
Qg(candidate_nodes,:) >= -0.04*repmat(x_es,1,T);
];

%%
%支路电流约束-加上这个约束运算速度过慢
Constraints = [Constraints, 0 <= Iij,Iij <= 1.1*Zij];
%支路功率约束
Constraints = [Constraints, -1.1*Zij <= Pij,Pij <= 1.1*Zij];
% 不能给支路无功功率施加约束
Constraints = [Constraints, -1.1*Zij <= Qij,Qij <= 1.1*Zij];
% 负荷损失的状态变量
Constraints = [Constraints,0<=lamda,lamda <= 1];
% 定义负荷重要程度
Importance = [10 1 3 1 1 1 10 3 1 1 1 3 1 3 1 1 1 10 1 1 3 1 1 1 10 1 1 1 3 1 10 1 1];
% 定义各种成本(单位:美元/千瓦时=万美元/10MWh)
pload_cost = 3.5;
ploss_cost = 0.08;
Pes_dc_cost = 0.2; % 储能放电成本 

%% 4.设目标函数-内层的目标函数-损耗成本+建设费用(单位:万美元)
% objective = sum(sum(Iij.*(r*ones(1,T)))) + sum(sum(pload))+sum(sum(-lamda.*pload))+;%网损最小+负荷损失最小;可以在sum(sum(-lamda.*pload))前乘以负荷重要性矩阵
ob1 = pload_cost * sum(sum((pload-lamda.*pload).*(Importance'*ones(1,4))));
ob2 = ploss_cost * sum(sum(Iij.*(r*ones(1,T))));
ob3 = Pes_dc_cost*sum(sum(Pes_dc));
% objective = ob1+ob2+ob3;

% 外层目标函数
Ees_cost = 105; % 电池配置价格 单位:万美元/10MWh 
w_ess = 0.003; % 储能成本权重-每年可能发生的 1-2 次 电网断电事故;储能预期使用年限：10 – 15 年
X = w_ess*Ees_cost*sum(Ees_max);

objective = ob1+ob2+ob3+X;
r_load = sum(sum(lamda.*pload))/sum(sum(pload)); % 负荷恢复率
V_bias = sum(abs(sqrt(U)-1))/33; % 电压偏差情况（分为四个周期）
R_load = sum((lamda.*pload)./pload,2)/4;

% 根据ChatGPT概念 R_load 节点负荷恢复率公式应为 
% R_load = sum((lamda.*pload),2)./sum(pload,2) % 因为在不同时刻下的节点负荷分母不同
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
% 目标函数值
ob1 = value(ob1);ob2 = value(ob2);ob3 = value(ob3);X = value(X); 
objective = value(objective);
r_load = value(r_load);
V_bias = value(V_bias);
R_load = value(R_load);
%% 6.输出AMPL模型
%saveampl(Constraints,objective,'mymodel');
%% 7.分析错误标志
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end

x_es_value = round(value(x_es));
x_es_T = value(x_es_T);

% 由CPLEX确定实际安装位置
es_nodes = candidate_nodes(x_es_value > 0.5);

% 保留所有候选节点容量
Ees_max_all = value(Ees_max);

% 仅输出已安装储能的容量
Ees_max = Ees_max_all(x_es_value > 0.5);

% 实际安装数量
N_es = sum(x_es_value);

%% 8.得到运行结果
U = value(U);Iij = value(Iij);Pij = value(Pij);Qij = value(Qij);Pg = value(Pg); Qg = value(Qg);
Zij = value(Zij);lamda = value(lamda);
% Fij = value(Fij);Wj = value(Wj);
Pes_c = value(Pes_c);Pes_dc = value(Pes_dc);Ees = value(Ees);
end