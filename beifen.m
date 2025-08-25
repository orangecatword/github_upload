%% 初始化变量
% 二阶锥约束的计算电压
U_cone = sdpvar(32,24,25,'full'); % 不同场景下的节点电压
% Gij = g; % 电导
% Bij = b; % 电纳
% theta_ij = zeros(32, 1);
% theta_ij = atan2(b, g); % 角度,用反正切函数计算
% theta_ij_expanded = repmat(reshape(theta_ij, [1, 32, 1]), [25, 1, 24]);

%% 支路链接关系
% 定义节点连接关系(示例,需根据实际拓扑调整)
% branches = {
%     33, 1;
%     1, 2;
%     2, 3;
%     3, 4;
%     4, 5;
%     5, 6;
%     6, 7;
%     7, 8;
%     8, 9;
%     9, 10;
%     10, 11;
%     11, 12;
%     12, 13;
%     13, 14;
%     14, 15;
%     15, 16;
%     16, 17;
%     1, 18
%     18, 19;
%     19, 20;
%     20, 21;
%     2, 22;
%     22, 23;
%     23, 24;
%     5, 25;
%     25, 26;
%     26, 27;
%     27, 28;
%     28, 29;
%     29, 30;
%     30, 31;
%     31, 32;
% };

%% 约束条件
    % 1.最最重要的风电和光伏出力约束
    % tic
    % % 维度扩展
    % LC_wt_expanded = repmat(LC_wt, [25, 1, 24]);
    % LC_pv_expanded = repmat(LC_pv, [25, 1, 24]);
    % cv = [cv;
    %     P_wt >= 0
    %     P_pv >= 0
    %     P_wt <= LC_wt_expanded
    %     P_pv <= LC_pv_expanded
    %     ];
    % toc


    % 3.功率平衡约束
    % 3.1 潮流约束
    % for hour = 1:24
    %     for s = 1:N_s
    %         for i = 1:N_n
    %     % 有无功功率平衡约束
    %         sum_P = 0;
    %         sum_Q = 0;
    %         for k = 1:N_n-1
    %             j = branches{k, 2};
    %             % if j == i
    %                 % sum_P = sum_P + U(s,branches{k,1}) * (Gij(k) * cos(theta_ij(k)) + Bij(k) * sin(theta_ij(k))); % 导纳值可以计算出;theta_ij(s)表示阻抗角
    %                 % sum_Q = sum_Q + U(s,branches{k,1}) * (Gij(k) * sin(theta_ij(k)) - Bij(k) * cos(theta_ij(k)));
    %             if branches{k,1} == i
    %                 sum_P = sum_P + U(s,branches{k,2},hour) * (Gij(k) * cos(theta_ij(k)) + Bij(k) * sin(theta_ij(k)));
    %                 sum_Q = sum_Q + U(s,branches{k,2},hour) * (Gij(k) * sin(theta_ij(k)) - Bij(k) * cos(theta_ij(k)));
    %             end
    %         end
    %         cv = [cv; U(s,i,hour) * sum_P == P(s,i,hour)];
    %         cv = [cv, U(s,i,hour) * sum_Q == Q(s,i,hour)];
    %         end
    %     end
    % end
    
    % 3.1 电压潮流约束(凸松弛)
    % 预计算所有支路的三角函数值
    % tic
    % Rij_expanded = repmat(reshape(r, 1, 32), [25, 1, 24]); % 电阻矩阵扩展到(25,32,24)
    % Xij_expanded = repmat(reshape(x, 1, 32), [25, 1, 24]); % 电抗矩阵扩展到(25,32,24)
    % 
    % for k = 1:N_n-1
    %     i = branches{k,1};
    %     j = branches{k,2};
    %     cv = [cv;
    %     U(:,j,:) == U(:,i,:) - 2*(Rij_expanded(:,k,:).*Pij(:,k,:) + Xij_expanded(:,k,:).*Qij(:,k,:)) + (Rij_expanded(:,k,:).^2 + Xij_expanded(:,k,:).^2).*Iij(:,k,:); 
    %     ];
    % end
    % toc

    % 3.2 有无功功率的潮流等式（凸松弛）(1.0版
    % 预计算节点连接关系
    % tic
    % node_in = cell(N_n-1,1);
    % node_out = cell(N_n-1,1);
    % branches = cell2mat(branches); 
    % for k = 1:32
    %     j = branches(k,2);
    %     % node_in{j} = find(branches(:,2) == j);  % 流入节点j的支路索引-流入的只可能有一个节点
    %     node_out{j} = find(branches(:,1) == j); % 流出节点j的支路索引
    %     if ~isempty(node_out{j})
    %     % cv = [cv;
    %     % sum(Qij(:,j,:) - Xij_expanded(:,j,:).*Iij(:,j,:), 2) == sum(Qij(:,node_out{j},:), 2) + Q(:,j,:);
    %     % sum(Pij(:,j,:) - Rij_expanded(:,j,:).*Iij(:,j,:), 2) == sum(Pij(:,node_out{j},:), 2) + P(:,j,:)];
    %     Q(:,j,:) = sum(Qij(:,j,:) - Xij_expanded(:,j,:).*Iij(:,j,:), 2) - sum(Qij(:,node_out{j},:), 2);
    %     P(:,j,:) = sum(Pij(:,j,:) - Rij_expanded(:,j,:).*Iij(:,j,:), 2) - sum(Pij(:,node_out{j},:), 2);
    %     else
    %     % sum(Qij(:,j,:) - Xij_expanded(:,j,:).*Iij(:,j,:), 2) == Q(:,j,:);
    %     % sum(Pij(:,j,:) - Rij_expanded(:,j,:).*Iij(:,j,:), 2) == P(:,j,:)];
    %     Q(:,j,:) = sum(Qij(:,j,:) - Xij_expanded(:,j,:).*Iij(:,j,:), 2)
    %     P(:,j,:) = sum(Pij(:,j,:) - Rij_expanded(:,j,:).*Iij(:,j,:), 2);
    %     end
    % end
    % toc
    % 
    % % 3.3 节点接入的各个功率等式-不知道Pdr的原理 暂时将这部分隐去(1.0版
    % tic
    % P_DG = P_wt + P_pv;
    % % 计算所有节点的净有功功率
    % for j = 1:N_n
    %     if j ~= 33 % 若非平衡节点
    %         cv = [cv; 
    %              P(:,j,:) == P_DG(:,j,:) - P_L(:,j,:) 
    %              ];
    %     else
    %         cv = [cv;
    %              P(:,33,:) == P_DG(:,33,:) - P_L(:,33,:) + P_en(:,1,:) 
    %              ];
    %     end
    % end
    % toc
    % 
    % % 3.4 二阶锥约束
    % tic
    % for k = 1:32  % 32条支路
    %     i = branches(k,1);  % 支路起始节点
    %     U_cone(:,k,:) = U(:,i,:);
    % end
    % Pij_vec = reshape(Pij, [], 1);
    % Qij_vec = reshape(Qij, [], 1);
    % Iij_vec = reshape(Iij, [], 1);
    % U_cone_vec = reshape(U_cone, [], 1);
    % % cv = [cv, cone([2*Pij_vec; 2*Qij_vec; Iij_vec-U_cone_vec], Iij_vec+U_cone_vec)];
    % cv = [cv; norm([2*Pij_vec; 2*Qij_vec; Iij_vec-U_cone_vec]) <= Iij_vec+U_cone_vec];
    % cv = [cv; cone([Pij_vec; Qij_vec], 7000)];
    % toc
 
    % 3.2 无功功率的潮流等式
    % for hour = 1:24
    %     for s = 1:N_s
    %         for j = 1:N_n
    %             % 计算流入节点j的净无功功率
    %             sum_Q_in = 0;
    %             sum_Q_out = 0;
    %             for k = 1:32
    %                 if branches{k,2} == j
    %                     sum_Q_in = sum_Q_in + (Qij(s,k,hour) - Xij(k)*Iij(s,k,hour)^2);
    %                 end
    %                 if branches{k,1} == j
    %                     sum_Q_out = sum_Q_out + Qij(s,k,hour);
    %                 end
    %             end
    %             cv = [cv; 
    %             sum_Q_out - sum_Q_in == Q(s,j,hour)
    %             ];
    %         end
    %     end
    % end

    % 3.3 有功功率的潮流等式
    % P_DG = P_wt + P_pv;
    % for hour = 1:24
    %     for s = 1:N_s
    %         for j = 1:N_n
    %             % 计算流入节点j的净有功功率
    %             sum_P_in = 0;
    %             sum_P_out = 0;
    % 
    %             % 计算∑[P_ij - R_ij*I_ij^2]
    %             for k = 1:32
    %                 if branches{k,2} == j
    %                     sum_P_in = sum_P_in + (Pij(s,k,hour) - Rij(k)*Iij(s,k,hour)^2);
    %                 end
    %                 if branches{k,1} == j
    %                     sum_P_out = sum_P_out + Pij(s,k,hour);
    %                 end
    %             end
    %             if j == 33 % 假设节点33是连接上级电网的节点
    %             % 节点j的有功功率平衡约束
    %             cv = [cv; 
    %                 sum_P_out - sum_P_in == P(s,j,hour);
    %                 P(s,j,hour) == P_DG(s,j,hour) - P_L(s,j,hour) + P_en(s,j,hour) - P_dr(s,j,hour)
    %             ];
    %             else
    %             cv = [cv; 
    %                 sum_P_out - sum_P_in == P(s,j,hour);
    %                 P(s,j,hour) == P_DG(s,j,hour) - P_L(s,j,hour) - P_dr(s,j,hour)
    %             ];
    %             end
    %         end
    %     end
    % end


    % 3.4 二阶锥约束(1.0版
    % for hour = 1:24
    %     for s = 1:25
    %         for k = 1:32  % 32条支路
    %             i = branches{k,1};  % 支路起始节点
    %             cv = [cv;
    %                  cone([2*Pij(s,k,hour); 2*Qij(s,k,hour); Iij(s,k,hour)-U_cone(s,k,hour)], Iij(s,k,hour)+U_cone(s,k,hour))];
    %             ];
    %         end
    %     end
    % end

    % 4.支路容量约束-已转化为支路有功与无功的二阶锥约束
    % Sij(j)和S_max(j) j分别表示线路j的视在功率及其上限(单位:MVA)-变成kVA后乘1000
    % 修改数据格式 或者放在最上面赋值即可
    % tic
    % Sij = Pij .* (1./(cos(theta_ij_expanded) * 1000));
    % cv = [cv; 
    % Sij <= 7000;
    % Sij >= 0;
    % ];  
    % toc

    % 5.DG运行约束- 未给出切削率上限
    % P_DG: 节点i在场景s下的有源出力上限
    % omega_DG: 分布式电源的切削率-切削率光伏出力除以光伏容量
    % P_DG_max节点总容量
    % omega_DG = P_DG .*(1./ LC_wt_expanded); % 当分母为0时赋值为0
    % 重点想法:在其他条件一致的情况下 光伏出力与光伏容量成正比关系
    % for i = 1:N_n
    %     for s = 1:N_s
    %         for hour = 1:24
    %         cv = [cv;
    %         % P_DG(s,i,hour) <= P_DG_max(s,i,hour) % 每个节点的DG装机容量上限为800kw
    %         P_DG <= 800;
    %         P_DG(s,i,hour) >= (1 - omega_DG(s,i,hour)) * 800  
    %         ];
    %         end
    %     end
    % end

    % 6.DG渗透率约束-见公式(21)-是上层的约束条件
    % 最大渗透率通常指的是配电网能够接纳的分布式电源（如光伏、风电等）的最大容量与电网最大负荷的比值
    % 最大渗透率为0.45
    % tic
    % P_DG_sum = sum(LC_wt+LC_pv);
    % 
    % cv = [cv;
    %     P_en <= (1-0.45)*P_DG_sum*(1/0.45);% 暂将发电量定义为不同场景下不同时间下 若乘以时间,就从功率变成功了
    %     ];
    % toc

%% 目标函数

% b.DG运行成本
% C_om = 0;
% for s = 1:N_s
%     for i = 1:N_n
%         for hour = 1:24
%         % C_om = C_om + t * (c_om_wt * sum(P_wt(s,i,:),3) + c_om_pv * sum(P_pv(s,i,:),3));
%         C_om = C_om + t(s) * (0.2 * sum(P_wt(s,i,:),3) + 0.2 * sum(P_pv(s,i,:),3)); % 0.2元/kwh
%         end
%     end
% end

% d.网络损失费用
% 公式表示对所有场景的网络损失费用进行求和,其中包含损失系数(c_loss)、损失功率(P_loss)和场景s对应天数(t)
% tic
% Rij_matrix = repmat(reshape(Rij, 1, 32), [25, 1, 24]); % 扩展电阻维度
% for s = 1:25  % 遍历所有场景
%     for k = 1:32
%         for hour = 1:24
%             P_loss(s) = P_loss(s)+Iij(s,k,hour)*Rij_matrix(s,k,hour);
%         end
%     end
% end
% 
% c_loss = 0.5;
% C_loss = sum(c_loss .* P_loss .* t);
% toc

%% 求解提速
% tic
% % 1. 基础求解器配置
% ops = sdpsettings('solver', 'cplex', 'verbose', 0);      
% 
% % 2. 计算资源优化
% ops.cplex.timelimit = 600;      % 延长至10分钟(原300秒)
% ops.cplex.parallel = 1;        % 启用并行计算(-1表示自动检测核心数)
% ops.cplex.threads = 0;         % 0表示自动分配线程
% 
% % 3. 求解精度控制  
% ops.cplex.mip.tolerances.mipgap = 0.03;  % 最优间隙从2%放宽到3%
% ops.cplex.mip.tolerances.absmipgap = 1; % 绝对间隙阈值(kW)
% ops.cplex.emphasis.numerical = 1;       % 增强数值稳定性
% 
% % 4. 搜索策略优化
% ops.cplex.mip.strategy.search = 2;      % 2=自动选择(平衡速度与质量)
% ops.cplex.mip.strategy.nodeselect = 3;  % 3=最优估计搜索
% ops.cplex.mip.strategy.heuristicfreq = 50; % 每50节点启发式搜索
% ops.cplex.mip.strategy.variableselect = 4; % 4=伪成本分支
% 
% % 5. 预求解增强
% ops.cplex.preprocessing.presolve = 1;    % 启用预求解器
% ops.cplex.preprocessing.reduce = 3;     % 最大程度简化模型
% ops.cplex.preprocessing.linear = 1;    % 线性化处理
% 
% % 6. 热启动配置(需历史解)
% if exist('prev_sol','var')
%     ops.cplex.advance = 2;              % 2=使用初始解
%     ops.cplex.start = prev_sol;         % 注入历史解
% end
% 
% % 执行优化
% optimize(cv, f, ops);
% toc