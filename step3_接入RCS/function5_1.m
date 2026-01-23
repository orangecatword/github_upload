% clear 
% clc
% warning off
%% 按照文章统一参数
function[objective,r_load,V_bias] = function5_1(RCS)
% function[objective,Psum_loss,Psum_load,U] = function5_1(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3)

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
for i=1:32
    upstream(i+1,i)=1;
end
upstream(21,33)=1;%支路33为21-8支路，流入节点21
upstream(15,34)=1;%支路34为15-9支路，流入节点15
upstream(22,35)=1;%支路35为22-12支路，流入节点22
upstream(33,36)=1;%支路36为33-18支路，流入节点33
upstream(29,37)=1;%支路37为29-25支路，流入节点29
for i=[1:17,19:21,23:24,26:32]
    dnstream(i,i)=1;
end
dnstream(2,18)=1;
dnstream(3,22)=1;
dnstream(6,25)=1;
% 5条流入，对应5条流出
dnstream(8,33)=1;
dnstream(9,34)=1;
dnstream(12,35)=1;
dnstream(18,36)=1;
dnstream(25,37)=1;
Umax=[1*1*ones(1,T);1.06*1.06*ones(32,T)];
Umin=[1*1*ones(1,T);0.94*0.94*ones(32,T)];

%% 1. 重新定义物理出力限值(Gemini新加入)
% 初始化所有节点出力上限为 0
Pgmax = zeros(nb, T); 
Qgmax = zeros(nb, T);

% (1) 主网节点 (节点 1)
Pgmax(1, :) = 1.0; % 对应 10MW
Qgmax(1, :) = 1.0; 

% (2) 定义分布式电源 (DER) 的物理容量
% 建议容量设得略大于风光出力的峰值，以防逆变器限功率
% 根据Case33bw,参照之前的分析：风机峰值约 1.5MW (0.15pu)，光伏峰值约 1MW (0.1pu)
pv_nodes = [7, 27];
wt_nodes = 12;
% es_nodes = [3, 22,32];
es_nodes = [18,24,30];
Pgmax(pv_nodes, :) = 0.2; 
Qgmax(pv_nodes, :) = 0.2; % 给予一定的无功补偿空间
Pgmax(wt_nodes, :) = 0.2; 
Qgmax(wt_nodes, :) = 0.2;
Pgmax(es_nodes, :) = 0.2; 
Qgmax(es_nodes, :) = 0.2;
%% 2.设变量
U = sdpvar(nb,T);%电压的平方
Iij = sdpvar(nl,T);%电流的平方
Pij = sdpvar(nl,T);%线路有功
Qij = sdpvar(nl,T);%线路无功
Pg = sdpvar(nb,T);%发电机有功
Qg = sdpvar(nb,T);%发电机无功
P = sdpvar(nb,T);
Q = sdpvar(nb,T);

% 定义节点状态变量-用以判断节点是否正常运行
u = binvar(nb,T);
lamda=sdpvar(33,T,'full');% 负荷损失的状态变量
Zij=binvar(nl,T);%网架结构               
% Z0=[ones(nl-nc,1);zeros(nc,1)];%初始拓扑                      
% assign(Zij,Z0);          

% 确认风光的安装位置 风: 12 光: 7 27
% Loc_pv_initial = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
% Loc_wt_initial = [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% Ppv = (Loc_pv_initial/1)' * (mpc.pv(9:12)/10); % 计算得到标幺值
% Pwt = (Loc_wt_initial/1.5)' * (mpc.wind(9:12)/10); % 除1.5是把风光机组拉到同一容量单位

% 储能状态
Ies_c = binvar(3, T,'full');% 充电状态
Ies_dc = binvar(3, T,'full');% 放电状态
% Ies_c1 = binvar(1, T,'full');% 充电状态
% Ies_dc1 = binvar(1, T,'full');% 放电状态
% Ies_c2 = binvar(1, T,'full');% 充电状态
% Ies_dc2 = binvar(1, T,'full');% 放电状态
% Ies_c3 = binvar(1, T,'full');% 充电状态
% Ies_dc3 = binvar(1, T,'full');% 放电状态

% 储能容量-根据论文修改
% Ees_max = Ees_max';
% Ees_max = [1.0188;      1.0663;     0.84204];
Ees_max = ones(3,1);
% Ees_max1 = 1; % 标幺值 实际为10MWh
% Ees_max2 = 1;
% Ees_max3 = 1;
% 储能设备剩余能量
Ees = sdpvar(3, T+1,'full'); 
% Ees1 = sdpvar(1, T+1,'full'); 
% Ees2 = sdpvar(1, T+1,'full');
% Ees3 = sdpvar(1, T+1,'full');

% 确认储能的安装位置  3  22 32
% Lc_pes1 = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% Lc_pes2 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
% Lc_pes3 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];

% 储能功率
Pes_c = sdpvar(3, T,'full');% 充电功率
Pes_dc = sdpvar(3, T,'full');% 放电功率
% Pes_c1 = sdpvar(1, T,'full');% 充电功率
% Pes_dc1 = sdpvar(1, T,'full');% 放电功率
% Pes_c2 = sdpvar(1, T,'full');% 充电功率
% Pes_dc2 = sdpvar(1, T,'full');% 放电功率
% Pes_c3 = sdpvar(1, T,'full');% 充电功率
% Pes_dc3 = sdpvar(1, T,'full');% 放电功率
%% 3.设约束
Constraints = [];    
%% 储能约束-可以取消充放电的功率和最好一致的约束，以增强储能供电能力
Constraints = [Constraints;
    Ies_c == 0;
    Ies_dc == 1; % 储能状态约束
    % sum(Pes_c) - sum(Pes_dc) == 0; % 充放电的功率和最好一致
    Pes_c == 0;
    Pes_c <= 0.1*Ies_c.*(Ees_max*ones(1,T)); % 储能容量的10% 
    Pes_dc >= 0;
    Pes_dc <= 0.1*Ies_dc.*(Ees_max*ones(1,T));
    Ees(:,1) == 0.5*Ees_max;
    Ees >= 0.2 * Ees_max*ones(1,T+1);
    Ees <= 0.8 * Ees_max*ones(1,T+1)
    ];
Constraints = [Constraints;
    Ees(:, 2:T+1) == Ees(:, 1:T) + 0.9*Pes_c - 1.1*Pes_dc];

% Constraints = [Constraints;
%     Ies_c1 + Ies_dc1 <= 1; % 储能状态约束
%     % sum(Pes_c1) - sum(Pes_dc1) == 0; % 充放电的功率和最好一致
%     Pes_c1 >= 0;
%     Pes_c1 <= Ies_c1 * 0.1*Ees_max1; % 储能容量的10% 
%     Pes_dc1 >= 0;
%     Pes_dc1 <= Ies_dc1 * 0.1*Ees_max1;
%     Ees1(1) == 0.5*Ees_max1;
%     Ees1 >= 0.2 * Ees_max1;
%     Ees1 <= 0.8 * Ees_max1
%     ];
% Constraints = [Constraints;
%     Ees1(2:T+1)==Ees1(1:T)+0.9*Pes_c1-1.1*Pes_dc1];

% Constraints = [Constraints;
%     Ies_c2 + Ies_dc2 <= 1; % 储能状态约束
%     % sum(Pes_c2) - sum(Pes_dc2) == 0; % 充放电的功率和最好一致
%     Pes_c2 >= 0;
%     Pes_c2 <= Ies_c2 * 0.1*Ees_max2; % 储能容量的10%
%     Pes_dc2 >= 0;
%     Pes_dc2 <= Ies_dc2 * 0.1*Ees_max2;
%     Ees2(1) == 0.5*Ees_max2;
%     Ees2 >= 0.2 * Ees_max2;
%     Ees2 <= 0.8 * Ees_max2
%     ];
% Constraints = [Constraints;
%     Ees2(2:T+1)==Ees2(1:T)+0.9*Pes_c2-1.1*Pes_dc2];

% Constraints = [Constraints;
%     Ies_c3 + Ies_dc3 <= 1; % 储能状态约束
%     % sum(Pes_c3) - sum(Pes_dc3) == 0; % 充放电的功率和最好一致
%     Pes_c3 >= 0;
%     Pes_c3 <= Ies_c3 * 0.1*Ees_max3; % 储能容量的10%
%     % Pes_c3<= 0;
%     Pes_dc3 >= 0;
%     Pes_dc3 <= Ies_dc3 * 0.1*Ees_max3;
%     Pes_dc3 >= 0.02;
%     Pes_dc3 <= 0.03;
%     Ees3(1) == 0.5*Ees_max3;
%     Ees3 >= 0.2 * Ees_max3;
%     Ees3 <= 0.8 * Ees_max3
%     ];
% Constraints = [Constraints;
%     Ees3(2:T+1)==Ees3(1:T)+0.9*Pes_c3-1.1*Pes_dc3];

%% 潮流约束
% 节点功率约束
Constraints = [Constraints, pload.*lamda - Pg== upstream*Pij - upstream*(Iij.*(r*ones(1,T))) - dnstream*Pij];
Constraints = [Constraints, qload.*lamda - Qg == upstream*Qij - upstream*(Iij.*(x*ones(1,T))) - dnstream*Qij];%节点注入无功
% Constraints = [Constraints, P == pload.*lamda - Pg+ Lc_pes1' * Pes_c1 - Lc_pes1' * Pes_dc1 + Lc_pes2' * Pes_c2 - Lc_pes2' * Pes_dc2 + Lc_pes3' * Pes_c3 - Lc_pes3' * Pes_dc3];
% Constraints = [Constraints, Q == qload.*lamda - Qg];

m = 1.06*1.06 - 0.94*0.94;
M = (ones(nl,T) - Zij)*m;             
Constraints = [Constraints, U(branch(:,1),:) - U(branch(:,2),:) <= M + 2*(r*ones(1,T)).*Pij + 2*(x*ones(1,T)).*Qij - ((r.^2 + x.^2))*ones(1,T).*Iij];
Constraints = [Constraints, U(branch(:,1),:) - U(branch(:,2),:) >= -M + 2*(r*ones(1,T)).*Pij + 2*(x*ones(1,T)).*Qij - ((r.^2 + x.^2))*ones(1,T).*Iij];
% 如果节点在线，电压在常规范围；如果节点离线，电压固定在 1.0。这能极大地辅助求解器收敛
V_nominal = 1.0^2;
Constraints = [Constraints, u*(0.94^2) + (1-u)*V_nominal <= U,U <= u*(1.06^2) + (1-u)*V_nominal];
%节点电压约束               
Constraints = [Constraints, U(1, :) == 1.0]; % 首节点电压为1   

%% 商品流约束-问题：出现了环流情况
% 论文观点为,除开故障状况,否则32条主干支路无法随便关断
pg_st = [1 7 12 27 es_nodes];
Fij=sdpvar(37,T,'full'); 
Wj=sdpvar(7,T,'full'); 
M_2=50;

A = 1:32;        % 定义集合A: [1,2,3,...,33]
% RCS = [14,20 21];
B = A(~ismember(A, RCS));  % 获取不属于RCS的元素
C = B(~ismember(B, [15,32]));  
% C = A(~ismember(A, [15,32]));  
for t=1:T
    for k=1:33
        if ~ismember(k,pg_st)
            node_out=branch(:,1)==k;
            node_in=branch(:,2)==k;
            Constraints=[Constraints,sum(Fij(node_out,t))-sum(Fij(node_in,t))==-u(k,t)]; 
            Constraints = [Constraints, u(k,t) >= Iij(node_out, t)];
            Constraints = [Constraints, u(k,t) >= Iij(node_in, t)];
        else
            node_out=branch(:,1)==k;
            node_in=branch(:,2)==k;
            Constraints=[Constraints,sum(Fij(node_out,t))-sum(Fij(node_in,t))==Wj(pg_st==k,t)-u(k,t)];
            Constraints = [Constraints, u(pg_st,t) == 1]; 
        end     
    end
        Constraints=[Constraints,Zij(15,t) == 0];
        Constraints=[Constraints,Zij(32,t) == 0]; 
        % Constraints=[Constraints,Zij(RCS,t) >= 0]; 
        % Constraints=[Constraints,Zij(RCS,t) <= 1];
        Constraints=[Constraints,Zij(C,t) == 1];
        %% 如果是正常运行W 与 F自然而然就是整数，问题就是如何在非强制的情况下避免环路？
        % Constraints=[Constraints,Zij(34,t) == 0];  
end
% Constraints=[Constraints,sum(Zij,1) <= sum(u,1)-1]; % 新增防止环流的点边约束
Constraints=[Constraints,sum(Zij,1) <= 32];
% Constraints=[Constraints,sum(Zij(1:32,:),1) >= 27];
Constraints=[Constraints,-M_2.*Zij<=Fij,Fij<=M_2.*Zij];
% Constraints=[Constraints,-M_2.*(2-Zij)<=Fij,Fij<=M_2.*(2-Zij)];

% Constraints=[Constraints,Wj>=1];
for i = 1:size(pg_st, 2)
    node_idx = pg_st(i);
    % 只要该 DER 节点在线，其供应量 Wj 就必须至少为 1（供应自己）
    Constraints = [Constraints, Wj(i, :) >= u(node_idx, :)];
end
Constraints = [Constraints, 0 <= Wj,Wj <= M_2 .* u(pg_st, :)];
for t = 1:T
    Constraints = [Constraints, sum(Wj(:,t)) == sum(u(:,t))];
end

%二阶锥约束  
for i = 1:37
    for t = 1:T
        Constraints = [Constraints;
            cone([2*Pij(i,t); 2*Qij(i,t);Iij(i,t) - U(branch(i,1),t)], Iij(i,t) + U(branch(i,1),t))]; 
    end
end

%% 通用约束
%发电机功率约束           
% Constraints = [Constraints, -Pgmax.*u <= Pg,Pg <= u .* Pgmax,-Qgmax.*u <= Qg,Qg <= u .* Qgmax];
Constraints = [Constraints, 0 <= Pg,Pg <= u .* Pgmax,-Qgmax.*u <= Qg,Qg <= u .* Qgmax];
%% Gemini
% Constraints = [Constraints, Pg(pv_nodes, :) == Ppv(pv_nodes, :).*u(pv_nodes, :) ];
% Constraints = [Constraints, Pg(wt_nodes, :) == Pwt(wt_nodes, :).*u(pv_nodes, :) ];
% Constraints = [Constraints, Pg(pv_nodes, :) <= Pgmax(pv_nodes, :).*u(pv_nodes, :) ];
% Constraints = [Constraints, Pg(wt_nodes, :) <= Pgmax(pv_nodes, :).*u(wt_nodes, :) ];
Constraints = [Constraints, Pg(es_nodes, :) == Pes_dc];

%%
%支路电流约束-加上这个约束运算速度过慢
Constraints = [Constraints, 0 <= Iij,Iij <= 1.1*Zij];
%支路功率约束
Constraints = [Constraints, -1.1*Zij <= Pij,Pij <= 1.1*Zij];
% 不能给支路无功功率施加约束
Constraints = [Constraints, -1.1*Zij <= Qij,Qij <= 1.1*Zij];
% 负荷损失的状态变量
Constraints=[Constraints,0<=lamda,lamda <= u];
% 定义负荷重要程度-参考论文:刘佳昕_极端灾害下有功-无功协同优化的两阶段配电网韧性提升策略
% Importance = [1 1 5 2 1 2 1 1 1 5 5 2 1 1 2 1 1 2 2 1 2 2 2 5 1 1 1 1 1 2 5 1 2];
% 同时考虑停电负荷和电网损耗的单位成本(单位:美元/千瓦时=万美元/10MWh)
pload_cost = 167;
ploss_cost = 3;
%% 4.设目标函数-内层的目标函数
% objective = sum(sum(Iij.*(r*ones(1,T)))) + sum(sum(pload))+sum(sum(-lamda.*pload));%网损最小+负荷损失最小;可以在sum(sum(-lamda.*pload))前乘以负荷重要性矩阵
objective = ploss_cost * sum(sum(Iij.*(r*ones(1,T)))) + pload_cost * sum(sum((pload-lamda.*pload)));
% Psum_loss = sum(sum(Iij.*(r*ones(1,T))));
r_load = sum(sum(lamda.*pload))/sum(sum(pload)); % 负荷恢复率
V_bias = sum(abs(sqrt(U)-1))/33; % 电压偏差情况（分为四个周期）
%% 5.设求解器
ops = sdpsettings('solver', 'cplex', 'verbose', 0);
ops.cplex.preprocessing.presolve = 1; % 启用预处理
ops.cplex.workmem = 4096; 
ops.cplex.mip.strategy.search = 1;  % 使用动态搜索
ops.cplex.mip.tolerances.mipgap = 0.05;  % 放宽最优间隙
ops.cplex.mip.tolerances.absmipgap = 1e-3; % 绝对间隙
ops.cplex.mip.strategy.heuristicfreq = 50; % 调整启发式频率
ops.cplex.threads = 1;        % ⭐ 关键 限制每个 CPLEX 求解器仅用 1 线程
ops.cplex.parallel = 0;       % ⭐ 关键 禁用 CPLEX 内部并行
ops.cplex.nodefileind = 2;    % 启用节点压缩磁盘文件
sol=optimize(Constraints,objective,ops);

% 目标函数值
objective = 100*value(objective);
% 节点状态值
u = value(u);
%% 6.输出AMPL模型
%saveampl(Constraints,objective,'mymodel');
%% 7.分析错误标志
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%% 8.得到运行结果
U = value(U);Iij = value(Iij);Pij = value(Pij);Qij = value(Qij);Pg = value(Pg);
Zij = value(Zij);lamda = value(lamda);Fij = value(Fij);Wj = value(Wj);
Pes_c = value(Pes_c);Pes_dc = value(Pes_dc);Ees = value(Ees);
% Pes_c1 = value(Pes_c1);Pes_dc1 = value(Pes_dc1);Ees1 = value(Ees1);
% Pes_c2 = value(Pes_c2);Pes_dc2 = value(Pes_dc2);Ees2 = value(Ees2);
% Pes_c3 = value(Pes_c3);Pes_dc3 = value(Pes_dc3);Ees3 = value(Ees3);

% figure(1)
% [XX,YY] = meshgrid(1:4,1:33);
% mesh(XX,YY,U);
% xlabel('时刻(h)');
% ylabel('节点序号');
% zlabel('电压幅值（pu）');
% title('24小时节点电压图');
% 
% figure(2)
% [XX,YY] =meshgrid(1:4,1:37);
% mesh(XX,YY,double(Iij));
% xlabel('时刻(h)');
% ylabel('节点序号');
% zlabel('线路电流（pu）');
% title('24小时线路电流标幺值图');    
% 
% figure(3)
% [XX,YY] =meshgrid(1:4,1:37);
% mesh(XX,YY,Pij);
% xlabel('时刻(h)');           
% ylabel('节点序号');
% zlabel('线路有功功率（pu）');
% title('24小时线路有功功率标幺值图');    
% 
% figure(4)
% [XX,YY] =meshgrid(1:4,1:37);
% mesh(XX,YY,Qij);
% xlabel('时刻(h)');             
% ylabel('节点序号');          
% zlabel('线路无功功率（pu）');      
% title('24小时线路无功功率标幺值图');       
 
% figure(5)
% [XX,YY] =meshgrid(1:4,1:37);
% mesh(XX,YY,Qij);
% xlabel('时刻(h)');             
% ylabel('节点序号');          
% zlabel('线路无功功率（pu）');      
% title('24小时线路无功功率标幺值图'); 
end