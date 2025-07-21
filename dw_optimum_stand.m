% 输出值都是计算目标函数需要的参数值
% 大多数变量都需要扩展维度 在时刻和场景下 包括支路功率 电压等
% 位置变量U P Q Pij Qij Sij 不仅要找到约束关系,还要明白它们之间如何计算-也包括出现在约束条件中的所有变量都得有个关系式囊括其中-还可能变得维度是 P_en
% 二阶锥约束-得到电流 支路功率与电压的关系
% 已有潮流约束-节点电压 节点功率的关系---结合二阶锥约束,应该能把这三者都解出来
% 损失功率计算-I^2*R
% 功率约束的两个经典公式-支路功率,支路电流与节点功率的关系;节点功率与 分布式电源出力 负荷出力 购电功率以及节点可转移功率(想法:应该有正有负)的关系


% 运算速度过慢 需要改变多层嵌套方式
function [t, U, P, Q, P_L, P_dr, f, P_en, P_loss] = dw_optimum_stand(LC_wt, LC_pv, P_wt, P_pv , r, x, g, b, p_l)

yalmip('clear') % 通过此方法,能加速求解过程

% 初始化变量
N_s = 25;    
N_n = 33;

% U0=1.0609; % 起始点电压(归一化后)

% 分时电价
d = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.45,0.56,0.56,0.45,0.45,0.4,0.56,0.4,0.45,0.45,0.45,0.4,0.4,0.2]; % 电价 单位元/1KWh

% 初始化变量
% 节点电压
U = sdpvar(25,33,24); % 不同场景下的节点电压

% 节点功率
P = sdpvar(25, 33, 24); 
Q = sdpvar(25, 33, 24);

% 支路功率
Pij = sdpvar(25, 32, 24);
Qij = sdpvar(25, 32, 24);
% 支路视在功率
Sij = sdpvar(25, 32, 24); 
% 支路电流
Iij = sdpvar(25, 32, 24);

Rij = r; % 电阻
Xij = x; % 电抗
Gij = g; % 电导
Bij = b; % 电纳
theta_ij = sdpvar(32, 1);
theta_ij = atan2(b, g); % 角度,用反正切函数计算


% 分布式电源的有功出力上限
P_DG = sdpvar(25, 33, 24);

% 分布式电源的切削率
omega_DG = sdpvar(25, 33, 24);

% 节点负荷功率
P_L = sdpvar(25, 33, 24);
% p_l
% p_l = repmat(p_l, 25, 1);
% p_l
% size(P_L)
% size(p_l)
% P_L(:,:,1) = p_l;

% 可转移负载功率 -是节点i在场景s下可转移负荷的总量
P_dr = sdpvar(25, 33, 24);

% 场景s在时刻t从上级电网购买电力的有功功率
P_en = sdpvar(25, 33, 24); 

% 损失功率-场景s下的总网络损耗功率
P_loss = sdpvar(1, 25);

% 定义节点连接关系(示例,需根据实际拓扑调整)
branches = {
    33, 1;
    1, 2;
    2, 3;
    3, 4;
    4, 5;
    5, 6;
    6, 7;
    7, 8;
    8, 9;
    9, 10;
    10, 11;
    11, 12;
    12, 13;
    13, 14;
    14, 15;
    15, 16;
    16, 17;
    1, 18
    18, 19;
    19, 20;
    20, 21;
    2, 22;
    22, 23;
    23, 24;
    5, 25;
    25, 26;
    26, 27;
    27, 28;
    28, 29;
    29, 30;
    30, 31;
    31, 32;
};


% 限制条件
    cv = [];

    % 维度扩展
    LC_wt_expanded = repmat(LC_wt, [25, 1, 24]);
    LC_pv_expanded = repmat(LC_pv, [25, 1, 24]);
    % 0.最最重要的风电和光伏出力约束
    cv = [cv;
        P_wt >= 0
        P_pv >= 0
        P_wt <= LC_wt_expanded
        P_pv <= LC_pv_expanded
        ];

    % 2.节点电压约束(标幺值 是已经除以12.66kv的结果)
    cv = [cv;
    % U(s,i) <= U_i_max
    % U(s,i) >= U_i_min
    % U <= 1.05
    % U >= 0.95
    % 若不除以标幺值
    U <= 1.05 * 12.66;  % 13.293 kV
    U >= 0.95 * 12.66;  % 12.027 kV
    ]; 

    % 3.功率平衡约束 - 有一个循环的 i Ns
    % 导纳值可以计算出
    % theta_ij(s)表示阻抗角

    % 3.1 潮流约束
    for hour = 1:24
        for s = 1:N_s
            for i = 1:N_n
        % 有无功功率平衡约束
            sum_P = 0;
            sum_Q = 0;
            for k = 1:N_n-1
                j = branches{k, 2};
                % if j == i
                    % sum_P = sum_P + U(s,branches{k,1}) * (Gij(k) * cos(theta_ij(k)) + Bij(k) * sin(theta_ij(k)));
                    % sum_Q = sum_Q + U(s,branches{k,1}) * (Gij(k) * sin(theta_ij(k)) - Bij(k) * cos(theta_ij(k)));
                if branches{k,1} == i
                    sum_P = sum_P + U(s,branches{k,2},hour) * (Gij(k) * cos(theta_ij(k)) + Bij(k) * sin(theta_ij(k)));
                    sum_Q = sum_Q + U(s,branches{k,2},hour) * (Gij(k) * sin(theta_ij(k)) - Bij(k) * cos(theta_ij(k)));
                end
            end
            cv = [cv; U(s,i,hour) * sum_P == P(s,i,hour)];
            cv = [cv, U(s,i,hour) * sum_Q == Q(s,i,hour)];
            end
        end
    end

    % 3.2 二阶锥约束
    % for hour = 1:24
    %     for s = 1:25
    %         for k = 1:32  % 32条支路
    %             i = branches{k,1};  % 支路起始节点
    %             cv = [cv;
    %                 Iij(s,k,hour)^2 == (Pij(s,k,hour)^2 + Qij(s,k,hour)^2) / U(s,i,hour)^2
    %             ];
    %         end
    %     end
    % end
    
    % 简化版二阶锥约束
    for hour = 1:24
        i_nodes = [branches{1:32,1}];  % 预计算所有支路起始节点
        cv = [cv, (Iij(:,:,hour).^2) == (Pij(:,:,hour).^2 + Qij(:,:,hour).^2) ./ (U(:,i_nodes,hour).^2)];
    end

    % 3.3 无功功率的潮流等式
    for hour = 1:24
        for s = 1:N_s
            for j = 1:N_n
                % 计算流入节点j的净无功功率
                sum_Q_in = 0;
                sum_Q_out = 0;
                
                for k = 1:32
                    if branches{k,2} == j
                        sum_Q_in = sum_Q_in + (Qij(s,k,hour) - Xij(k)*Iij(s,k,hour)^2);
                    end
                    if branches{k,1} == j
                        sum_Q_out = sum_Q_out + Qij(s,k,hour);
                    end
                end
                cv = [cv; 
                sum_Q_out - sum_Q_in == Q(s,j,hour)
                ];
            end
        end
    end

    % 3.4 有功功率的潮流等式
    P_DG = P_wt + P_pv;
    for hour = 1:24
        for s = 1:N_s
            for j = 1:N_n
                % 计算流入节点j的净有功功率
                sum_P_in = 0;
                sum_P_out = 0;
                
                % 计算∑[P_ij - R_ij*I_ij^2]
                for k = 1:32
                    if branches{k,2} == j
                        sum_P_in = sum_P_in + (Pij(s,k,hour) - Rij(k)*Iij(s,k,hour)^2);
                    end
                    if branches{k,1} == j
                        sum_P_out = sum_P_out + Pij(s,k,hour);
                    end
                end
                if j == 33 % 假设节点33是连接上级电网的节点
                % 节点j的有功功率平衡约束
                cv = [cv; 
                    sum_P_out - sum_P_in == P(s,j,hour);
                    P(s,j,hour) == P_DG(s,j,hour) - P_L(s,j,hour) + P_en(s,j,hour) - P_dr(s,j,hour)
                ];
                else
                cv = [cv; 
                    sum_P_out - sum_P_in == P(s,j,hour);
                    P(s,j,hour) == P_DG(s,j,hour) - P_L(s,j,hour) - P_dr(s,j,hour)
                ];
                end
            end
        end
    end

    % 4.支路容量约束
    % Sij(j)和S_max(j) j分别表示线路j的视在功率及其上限(单位:MVA)-变成kVA后乘1000
    % 修改数据格式 或者放在最上面赋值即可
    theta_ij_expanded = repmat(reshape(theta_ij, [1, 32, 1]), [25, 1, 24]);
    Sij = Pij./(cos((theta_ij_expanded)) * 1000);
    cv = [cv; 
    % Sij(j) <= S_max(j)
    Sij <= 7000;
    Sij >= 0;

    
    ];  

    % 6.DG运行约束
    % P_DG: 节点i在场景s下的有源出力上限
    % omega_DG: 分布式电源的切削率-切削率光伏出力除以光伏容量
    % P_DG_max节点总容量
    omega_DG = P_DG ./ LC_wt_expanded;
    for i = 1:N_n
        for s = 1:N_s
            for hour = 1:24
            omega_DG(s,i,hour) = P_DG(s,i,hour)/ LC_wt_expanded(s,i,hour);
            cv = [cv;
            % P_DG(s,i,hour) <= P_DG_max(s,i,hour) % 每个节点的DG装机容量上限为800kw
            P_DG <= 800;
            P_DG(s,i,hour) >= (1 - omega_DG(s,i,hour)) * 800  
            ];
            end
        end
    end
        cv = [
            cv;
            % omega_DG(i) <= omega_DG_max(i), % 可再生能源最高渗透率为45%
            omega_DG <= 0.45;
            omega_DG >= 0
            ];    


    % 7.可转移负荷量约束
    for i = 1:N_n
        for s = 1:N_s
            cv = [cv;
            P_dr(s,i) >= 0;          % P_dr(s,i) 是节点i在场景s下可转移负荷的总量
            % P_dr(s,i) <= P_dr_max(s,i)
            P_dr <= 37.15
            ]; 
        end
    end

    % 天数约束-由pk转为t
    % pk*365 再四舍五入，这里省去约束过程,即不需要pk参数了
    t = [10 10 10 11 15 25 25 25 27 38 10 10 10 11 16 9 9 9 9 13 11 12 11 12 17]; % 不同场景的运行天数


% 下层成本函数
% b.DG运行成本
% 计算总运行维护成本C_om
% N_s: 场景数量
% N_n: 节点数量
% t: 场景s的天数
% c_om_wt: 风电单位运维成本
% c_om_pv: 光伏单位运维成本
% P_wt: N_n*N_s矩阵,节点i场景s的风电功率
% P_pv: N_n*N_s矩阵,节点i场景s的光伏功率

C_om = 0;
for s = 1:N_s
    for i = 1:N_n
        for hour = 1:24
        % C_om = C_om + t * (c_om_wt * sum(P_wt(s,i,:),3) + c_om_pv * sum(P_pv(s,i,:),3));
        C_om = C_om + t(s) * (0.2 * sum(P_wt(s,i,:),3) + 0.2 * sum(P_pv(s,i,:),3)); % 0.2元/kwh
        end
    end
end

% c.上级电网购电成本
% 计算总环境成本C_en
% N_s: 场景数量
% d: 25*24矩阵,是系统的实际价格
% P_en: N_s*24矩阵,场景s在时刻t从上级电网购买电力的有功功率
% 看具体潮流约束公式,如果涉及到了单个节点约束,那么还是得扩展到三维
% t: 场景概率向量(1*N_s)
C_en = 0;
P_en_sum = sum(P_en, 2); 
for s = 1:N_s
    for hour = 1:24
        C_en = C_en + d(1,hour) .* P_en_sum(s,hour) .* t(s);
    end
end

% C_en = sum(sum(d .* P_en_sum .* t, 2));


% d.网络损失费用
% 公式表示对所有场景的网络损失费用进行求和,其中包含损失系数(c_loss)、损失功率(P_loss)和场景s对应天数(t)
for s = 1:25  % 遍历所有场景
    P_loss(s) = 0;  % 初始化场景s的损失功率
    for hour = 1:24  % 遍历所有小时
        for k = 1:32  % 遍历所有支路
            P_loss(s) = P_loss(s) + Iij(s,k,hour)^2 * Rij(k);
        end
    end
end
c_loss = 0.5;
C_loss = sum(c_loss .* P_loss .* t);

% e.实施需求响应的补偿成本
% P_dr(j,s)表示场景s下负载节点j处的可转移负载功率-可转移负载功率:配电网中能够通过需求响应(Demand Response, DR)机制在时间或空间上调整其用电行为的负荷功率
% c_dr表示每单位功率的可转移负载的补偿成本
C_dr = 0;
c_dr = 0.2;
P_dr_sum = sum(P_dr, 3); 
for s = 1:N_s
    for j = 1:N_n
        C_dr = C_dr + c_dr * P_dr_sum(s,j) * t(s);
    end
end


% 下层模型目标函数_DG总运营成本-除开投资成本外的年综合成本,之后再修改
f = C_om + C_en + C_loss + C_dr;

% 求解
ops = sdpsettings('solver', 'cplex', 'verbose', 0);
ops.cplex.timelimit = 300;  % 减少求解时间限制
ops.cplex.mip.tolerances.mipgap = 0.02;  % 放宽最优间隙
ops.cplex.parallel = 1;  % 启用并行计算
ops.cplex.mip.strategy.search = 1;  % 使用动态搜索
ops.cplex.mip.strategy.heuristicfreq = 100; % 调整启发式频率


optimize(cv, f, ops);

% 提取结果
t = value(t);
P_pv = value(P_pv);
P_wt = value(P_wt);
P_dr = value(P_dr);
f = value(f);
