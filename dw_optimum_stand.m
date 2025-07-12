% 输出值都是计算目标函数需要的参数值
function [P_wt, P_pv, P_DG, P_DG_max, U, P, Q, P_dr, f] = dw_optimum_stand(LC_wt, LC_pv,C_res, pk, g, b)


yalmip('clear') % 通过此方法,能加速求解过程


% 初始化变量
N_s = 25;    
N_n = 33;

U0=1.0609; % 起始点电压(归一化后)


% 初始化变量
% 电压约束
U = sdpvar(25,33); % 不同场景下的节点电压

P = sdpvar(25, 33); % 场景s时流入节点i的有功功率和无功功率
Q = sdpvar(25, 33);


% 导纳
Gij = sdpvar(32, 1); % 电导
Bij = sdpvar(32, 1); % 电纳
theta_ij = sdpvar(32, 1); % 角度

% 支路有无功功率
Pij = sdpvar(32, 1);
Qij = sdpvar(32, 1);
% 支路视在功率
Sij = sdpvar(32, 1); 

Gij = g;
Bij = b;
theta_ij = atan2(b, g) % 反正切函数
Sij = Pij./(cos((theta_ij)) * 1000) % 求视在功率


% 分布式电源的有功出力上限(某个场景的某个节点)
P_DG = sdpvar(25, 33);
P_DG_max = sdpvar(25, 33);
% 分布式电源的切削率
omega_DG = sdpvar(33, 1);
omega_DG_max = sdpvar(33, 1);

% 可转移负载功率
P_dr = sdpvar(25, 33);

% 不同场景的运行天数
t = sdpvar(1, 25);

% 节点的风电与光伏出力(某个场景的某个节点)
P_wt = sdpvar(25, 33); 
P_pv = sdpvar(25, 33);

% 定义节点连接关系(示例,需根据实际拓扑调整)
branches = {
    % 1, 2;
    % 2, [3,19];
    % 3, [4,23];
    % 4, 5;
    % 5, 6;
    % 6, [7,26];
    % 7, 8;
    % 8, 9;
    % 9, 10;
    % 10, 11;
    % 11, 12;
    % 12, 13;
    % 13, 14;       
    % 14, 15;       
    % 15, 16;
    % 16, 17;
    % 17, 18;       
    % 18, [];
    % 19, 20;
    % 20, 21;
    % 21, 22;        
    % 22, [];        
    % 23, 24;    
    % 24, 25;    
    % 25, [];    
    % 26, 27;         
    % 27, 28;         
    % 28, 29;         
    % 29, 30;    
    % 30, 31;         
    % 31, 32;          
    % 32, 33;
    % 33, [];
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

size(P_wt)
size(LC_wt)
% 限制条件
    cv = [];

    % 0.忘了最最重要的风电和光伏出力约束
    cv = [cv;
        P_wt >= 0
        P_pv >= 0
        P_wt <= repmat(LC_wt, 25, 1)
        P_pv <= repmat(LC_pv, 25, 1)
        ];

    % 2.节点电压约束(标幺值 是已经除以12.66kv的结果)
    cv = [cv,
    % U(s,i) <= U_i_max
    % U(s,i) >= U_i_min
    U <= 1.05
    U >= 0.95
    ]; 

    % 3.功率平衡约束 - 有一个循环的 i Ns
    % 导纳值可以计算出
    % P(s,i) Q(s,i)表示场景s时流入节点i的有功功率和无功功率
    
    % Gij Bij表示节点i和节点j之间导纳的真实的实部和虚部
    % theta_ij(s)表示阻抗角

    % 再修改
    for s = 1:numScenarios
    for i = 1:numNodes
        % 有功功率平衡约束
        sum_P = 0;
        for k = 1:numLines
            j = branches(k, 2);
            if j == i
                sum_P = sum_P + U(s,branches(k,1)) * (Gij(k) * cos(theta_ij(k)) + Bij(k) * sin(theta_ij(k)));
            elseif branches(k,1) == i
                sum_P = sum_P + U(s,branches(k,2)) * (Gij(k) * cos(theta_ij(k)) + Bij(k) * sin(theta_ij(k)));
            end
        end
        Constraints = [Constraints, U(s,i) * sum_P == P(s,i)];
        
        % 无功功率平衡约束
        sum_Q = 0;
        for k = 1:numLines
            j = branches(k, 2);
            if j == i
                sum_Q = sum_Q + U(s,branches(k,1)) * (Gij(k) * sin(theta_ij(k)) - Bij(k) * cos(theta_ij(k)));
            elseif branches(k,1) == i
                sum_Q = sum_Q + U(s,branches(k,2)) * (Gij(k) * sin(theta_ij(k)) - Bij(k) * cos(theta_ij(k)));
            end
        end
        Constraints = [Constraints, U(s,i) * sum_Q == Q(s,i)];
    end
end
    for s = 1:25
         for i = 1:33
             for j = 1:33 
    %             next = branches{node,1};
    %             current = branches{node,2};
    %     cv = [cv,  
    %             U(s,i) * sum(U(s,:).* (Gij(current).* cos(theta_ij(current)) + Bij(current).* sin(theta_ij(current)))) == P(s,i) % 有功功率平衡约束
    %             U(s,i) * sum(U(s,:).* (Gij(current).* sin(theta_ij(current)) - Bij(current).* cos(theta_ij(current)))) == P(s,i) % 无功功率平衡约束
    %     ];  
    %         end
    %     end
    % end


    % 4.支路容量约束
    % Sij(j)和S_max(j) j分别表示线路j的视在功率及其上限(单位:MVA)
    % Sij = Pij./(cos((theta_ij)) * 1000)
    cv = [cv; 
    % Sij(j) <= S_max(j)
    % 应该是已经归一化的结果
    Sij <= 7
    Sij >= 0
    ];  

    % 6.DG运行约束
    % P_DG: 节点i在场景s下的有源出力上限
    % omega_DG(i): 分布式电源的切削率
    for i = 1:N_n
        for s = 1:N_s
            cv = [cv,
            P_DG(s,i) <= P_DG_max(s,i)
            P_DG(s,i) >= (1 - omega_DG(i)) * P_DG_max(s,i)  
            ];
        end
        cv = [omega_DG(i) <= omega_DG_max(i),
            omega_DG(i) >= 0
            ];    
    end

    % 7.可转移负荷量约束
    for i = 1:N_n
        for s = 1:N_s
            cv = [cv, 
            0 <= P_dr(s,i),          % P_dr(s,i) 是节点i在场景s下的需求响应功率
            % P_dr(s,i) <= P_dr_max(s,i)
            P_dr(s,i) <= 37.15
            ]; 
        end
    end

    % 天数约束-由pk转为t
    % pk*365 再四舍五入，这里省去约束过程
    % cv = [cv, sum(t) == 365];
    % for i = 1:N_s
    %     cv = [cv, t(i) >= floor(365*pk(i))-1];
    %     cv = [cv, t(i) <= ceil(365*pk(i))+1];
    % end
    t = [10 10 10 11 15 25 25 25 27 38 10 10 10 11 16 9 9 9 9 13 11 12 11 12 17]; 


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
        % C_om = C_om + t * (c_om_wt * P_wt(s,i) + c_om_pv * P_pv(s,i));
        % C_om = C_om + t * (0.2 * P_wt(s,i) + 0.2 * P_pv(s,i)); 0.2元/kwh
        C_om = C_om + t(s) * 24 * (0.2 * P_wt(s,i) + 0.2 * P_pv(s,i)); % kw 和 kwh 单位换算 调整计算单位 一天24h
    end
end


% 下层模型目标函数_DG总运营成本-除开投资成本外的年综合成本,之后再修改
f = C_om;

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
