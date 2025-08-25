% 2025.8.21后一些思路:
% 1.简化了节点的可转移负荷功率;
% 4.参考熊壮壮的主程序约束-若按熊壮壮程序改 定义的变量大小应该 25 33 24变为33 25 24;以及定义upstream;

% 当在上层模型中完成选址定容后,下层模型的变量还剩下 可转移负荷出力以及上级电网给首端节点的功率
function [t, U, P, Q, P_dr_out, f, P_en, P_loss] = dw_optimum_stand(LC_wt, LC_pv, C_res, r, x, g, b, pload, qload,branch,pv,wt)

size(LC_wt) % 1    33
size(C_res) % 25    48
size(r) % 32     1
size(pload) % 33    24
% 目前不为0的变量:输出的变量,P_wt,P_pv,P_DG,theta_ij,t,d;把剩下的变量和约束都列出来 看看是不是约束少了
% 目前为0的变量:U, P, Q(Q第33列有非法值), P_L, P_dr, f, P_en, P_loss,Pij,Qij,Sij,Iij
% 至2025.8.19 不为0的变量增加P_L, P_dr,P_en
yalmip('clear') % 通过此方法,能加速求解过程

% 初始化变量
N_s = 25;    
N_n = 33;

% U0=1.0609; % 起始点电压(归一化后)

% 分时电价
d = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.45,0.56,0.56,0.45,0.45,0.4,0.56,0.4,0.45,0.45,0.45,0.4,0.4,0.2]; % 电价 单位元/1KWh

% 天数约束-由pk转为t; pk*365 再四舍五入，这里省去约束过程,即不需要pk参数了
t = [10 10 10 11 15 25 25 25 27 38 10 10 10 11 16 9 9 9 9 13 11 12 11 12 17]; % 不同场景的运行天数

% 初始化变量
% 节点电压的平方
U = sdpvar(33,24,25,'full'); % 不同场景下的节点电压

% 节点功率
P = sdpvar(33, 24, 25, 'full'); 
Q = sdpvar(33, 24, 25, 'full');

% 支路功率
Pij = sdpvar(32, 24, 25, 'full');
Qij = sdpvar(32, 24, 25, 'full');
% 支路视在功率
Sij = sdpvar(32, 24, 25, 'full'); 
% 支路电流的平方
Iij = sdpvar(32, 24, 25, 'full'); 

Rij = r; % 电阻
Xij = x; % 电抗

% 光伏和风电出力
P_wt = sdpvar(33, 24, 25, 'full');
P_pv = sdpvar(33, 24, 25, 'full');
% 分布式电源的有功出力
P_DG = sdpvar(33, 24, 25, 'full');
% 光伏风电出力与其容量成比
ratio_pv = zeros(33, 1);
ratio_wt = zeros(33, 1);

% 分布式电源的切削率
omega_DG = sdpvar(33, 24, 25, 'full');

% 节点负荷功率
P_L = repmat(reshape(pload, [33,24,1]), [1,1,25]);
Q_L = repmat(reshape(qload, [33,24,1]), [1,1,25]);

% 可转移负载功率 -是节点i在场景s下可转移负荷的总量
P_dr_out = sdpvar(33, 24, 25, 'full'); % 转移出的功率
P_dr_in = sdpvar(33, 24, 25, 'full'); % 转移入的功率
P_dr_sum = sdpvar(33, 25, 'full'); % 24时刻中正向可转移负荷之和
% 可转移负荷状态
u_out = binvar(33, 24, 25,'full');% 高峰转出
u_in = binvar(33, 24, 25,'full');% 低峰转入

% 场景s在时刻t从上级电网购买电力的有功功率-其实只是上级电网传到了33节点,大小应改变
P_en = sdpvar(1, 24, 25, 'full'); 

% 损失功率-场景s下的总网络损耗功率
P_loss = sdpvar(1, 25,'full');
P_loss_big = sdpvar(32, 24, 25, 'full');

% 目标函数
f = sdpvar(1, 1);

% 熊壮壮
% 网络拓扑构建
nb = 33;%节点数
nl = 32;%支路数
upstream = zeros(nb,nl);
dnstream = zeros(nb,nl);

% 为什么是从1到32 i到i本身值为1呢
for i = 1:nl
    upstream(i,i) = 1;
end

for i=[1:16,18:20,22:23,25:31] % 即存在子节点的节点集
    dnstream(i,i+1) = 1;
end

% 分支支路
dnstream(1,18) = 1;
dnstream(2,22) = 1;
dnstream(5,25) = 1;
dnstream(33,1) = 1;

% 完成节点风电光伏的赋值
tic
for i = 1:33
    if i == 13 || i == 17 || i == 25
        ratio_wt(i) = LC_wt(i)/1000; % LC_wt(i)/sum(LC_wt)*sum(LC_wt)/1000-计算风电节点的比例乘以容量的基准值
        % ratio_wt(i) = 1;
    elseif i == 4 || i == 7 || i == 27
        ratio_pv(i) = LC_pv(i)/sum(LC_pv); % 计算光伏节点的比例
        % ratio_pv(i) = 1;
    end
end

% 风电功率分配
% size(reshape(repmat(C_res(:,25:48)', [33,1]), [33,24,25]))
% size(repmat(reshape(ratio_wt,33,1,1), [1,24,25]))
% P_wt = 3*reshape(repmat(C_res(:,25:48)', [33,1]), [33,24,25]).*repmat(reshape(ratio_wt,33,1,1), [1,24,25]);

% 光伏功率分配
% P_pv = 3*reshape(repmat(C_res(:,1:24), [1,33]), [25,33,24]) .* repmat(reshape(ratio_pv,1,33,1), [25,1,24]);
% P_pv = reshape(repmat(C_res(:,1:24)', [33,1]), [33,24,25]).*repmat(reshape(ratio_pv,33,1,1), [1,24,25])
% 
% P_pv = repmat(reshape(pv, [1,24,1]), [33,1,25]).*repmat(reshape(ratio_pv,33,1,1), [1,24,25]);
% P_wt = repmat(reshape(wt, [1,24,1]), [33,1,25]).*repmat(reshape(ratio_wt,33,1,1), [1,24,25]);
size(reshape(repmat(C_res(:,1:24)', [33,1]), [33,24,25]))
size(repmat(reshape(ratio_pv,33,1,1), [1,24,25]))
P_pv = reshape(repmat(C_res(:,1:24)', [33,1]), [33,24,25]).*repmat(reshape(ratio_pv,33,1,1), [1,24,25]);
P_wt = reshape(repmat(C_res(:,25:48)', [33,1]), [33,24,25]).*repmat(reshape(ratio_wt,33,1,1), [1,24,25]);
toc

% 限制条件
    cv = [];

    U(33,:,:) = 1.03 ^2; % 首端节点电压标幺值
    % 2.节点电压约束(标幺值 是已经除以12.66kv的结果)-修改U的限制值后 变量无解 所以还是约束公式取值的问题
    cv = [cv;
    U <= 1.1025; %1.05平方1.1025  
    U >= 0.9025 %0.95平方0.9025
    ];

    cv = [cv, 0 <= Iij,Iij <= 6];%基于熊壮壮程序给电流一个约束

    % 7.可转移负荷量约束-可转移负荷是在不同时刻的负荷调动-节点t时刻可转移负荷量为正值,代表节点要在t时刻的原本负荷量减去可转移负荷量
    % 可转移负荷状态约束
    cv = [cv;
    u_out + u_in <=1 % 可转移负荷状态约束
    sum(P_dr_out,2) - sum(P_dr_in,2) == 0;
    P_dr_out >= 0;
    P_dr_out <= u_out*37.15;
    P_dr_in >= 0;
    P_dr_in <= u_in*37.15
    ];
    
    %% 3.熊壮壮版的潮流约束
    %% 潮流约束-判断公式具体是否正确
    % 节点净负荷计算公式
    % 省去了拓扑结构,仅用两行公式总结
    for j=1:25
        P(:,:,j) = upstream*Pij(:,:,j) - upstream*(Iij(:,:,j).*(r*ones(1,24))) - dnstream*Pij(:,:,j);
        Q(:,:,j) = upstream*Qij(:,:,j) - upstream*(Iij(:,:,j).*(x*ones(1,24))) - dnstream*Qij(:,:,j);
    end
    cv = [cv, Q == Q_L]; % 和大多数论文的约束公式对比,等式右侧整体加了个负号
    for j = 1:33
            if j ~= 33 % 若非平衡节点
                cv = [cv, P(j,:,:) == -(P_L(j,:,:) + P_dr_in(j,:,:) - P_dr_out(j,:,:) - P_pv(j,:,:)- P_wt(j,:,:))]; % 和大多数论文的约束公式对比,等式右侧整体加了个负号
            else
                cv = [cv, P(j,:,:) == -(P_L(j,:,:) + P_dr_in(j,:,:) - P_dr_out(j,:,:) - P_pv(j,:,:)- P_wt(j,:,:) - P_en(1,:,:))]; % 和大多数论文的约束公式对比,等式右侧整体加了个负号
            end
    end
    % 电压约束公式
    tic
    for j =1:25
        cv = [cv, U(branch(:,2),:,j) == U(branch(:,1),:,j) - 2*(r*ones(1,24)).*Pij(:,:,j) - 2*(x*ones(1,24)).*Qij(:,:,j) + ((r.^2+x.^2)*ones(1,24)).*Iij(:,:,j)];
        % cv = [cv, U(branch(:,1),:,j).*Iij(:,:,j) >= Pij(:,:,j).^2 + Qij(:,:,j).^2]; % 对熊壮壮论文(4-6)进行锥松弛
        % cv = [cv, 7 <= Pij(:,:,j).^2 + Qij(:,:,j).^2]; % 视在功率约束 7为标幺值
        Pij_vec = reshape(Pij, [], 1);
        Qij_vec = reshape(Qij, [], 1);
        Iij_vec = reshape(Iij, [], 1);
        U_vec = reshape(U(branch(:,2),:,:), [], 1);
        cv = [cv; norm([2*Pij_vec; 2*Qij_vec; Iij_vec-U_vec]) <= Iij_vec+U_vec];
        cv = [cv; cone([Pij_vec; Qij_vec], 7000)];
    end
    toc

% 下层成本函数
% b.DG运行成本
% 计算总运行维护成本C_om
% N_s: 场景数量
% N_n: 节点数量
% t: 场景s的天数
% c_om_wt: 风电单位运维成本
% c_om_pv: 光伏单位运维成本

% b.DG运行成本(简化版)
tic
% 32 24 25
% P_wt_sum = sum(P_wt, 2);  
% P_pv_sum = sum(P_pv, 2);
% C_om = C_om + sum(t.* (0.2 * sum(P_wt_sum, 1) + 0.2 * sum(P_pv_sum, 1)))  * 24;
% 向量化计算
C_om = sum(sum(0.2*P_wt,2)+sum(0.2*P_pv,2),1);% 大小为 1 1 25
C_om = C_om.*reshape(repmat(t', [1,1]), [1,1,25]);% 此时大小为1 1 25
% C_om = sum(sum(0.2*P_wt,2)+sum(0.2*P_pv,2),1).*t;
size(C_om)
C_om = sum(sum(C_om));
toc

% c.上级电网购电成本
% 计算总环境成本C_en
% N_s: 场景数量
% P_en: N_s*24矩阵,场景s在时刻t从上级电网购买电力的有功功率
% 看具体潮流约束公式,如果涉及到了单个节点约束,那么还是得扩展到三维
% t: 场景概率向量(1*N_s)
tic
size(P_en)
% P_en_sum = sum(P_en, 1); 
C_en = 0;
% C_en = C_en + d*sum(P_en, 1).*reshape(repmat(t', [1,1]), [1,1,25]);
for hour = 1:24
    for s = 1:25
        C_en = C_en + sum(d(1,hour) * P_en(1,hour,s) .* t(s));
    end
end
toc


% d.网络损失费用(简化版)
% 预计算电阻矩阵
tic
% 向量化计算
% P_loss = squeeze(sum(Iij .* Rij_expanded, 2));
for j = 1:25
    P_loss_big(:,:,j) = Iij(:,:,j).*(r*ones(1,24));
end
size(P_loss_big)
size(sum(sum(P_loss_big, 1),2))
P_loss = reshape(sum(sum(P_loss_big, 1),2), 1, 25);
size(P_loss)
   
c_loss = 0.5;
C_loss = sum(c_loss .* (t .* P_loss));
size(C_loss)
toc

% e.实施需求响应的补偿成本
% P_dr(j,s,hour)表示场景s下负载节点j处24小时的可转移负载功率-可转移负载功率:配电网中能够通过需求响应(Demand Response, DR)机制在时间或空间上调整其用电行为的负荷功率
% c_dr表示每单位功率的可转移负载的补偿成本

tic
% c_dr = 0.2;
% C_dr = sum(c_dr * P_dr_sum * t);

C_dr = sum(sum(0.2*P_dr_out,2),1);% 大小为 1 1 25
C_dr = C_dr.*reshape(repmat(t', [1,1]), [1,1,25]);% 此时大小为1 1 25
C_dr = sum(sum(C_dr));
size(C_dr)
toc

% 下层模型目标函数_DG总运营成本-除开投资成本外的年综合成本,之后再修改
size(C_om)
size(C_en)
size(C_loss)
size(C_dr)
f = C_om + C_en + C_loss + C_dr; % 不同项之间不能差太多
size(f)

% 求解
tic
ops = sdpsettings('solver', 'cplex', 'verbose', 0);
ops.cplex.timelimit = 300;  % 减少求解时间限制
ops.cplex.mip.tolerances.mipgap = 0.02;  % 放宽最优间隙
ops.cplex.parallel = 1;  % 启用并行计算
ops.cplex.mip.strategy.search = 1;  % 使用动态搜索
ops.cplex.mip.strategy.heuristicfreq = 100; % 调整启发式频率
% ops=sdpsettings('solver', 'cplex');
optimize(cv, f, ops);
toc

% 提取结果
t = value(t);
U = value(U);
Iij = value(Iij);
P = value(P);
Q = value(Q);
Pij = value(Pij);
Qij = value(Qij);
P_pv = value(P_pv);
P_wt = value(P_wt);
P_dr_out = value(P_dr_out);
size(f)
f = value(f);
C_en = value(C_en);
C_loss = value(C_loss);
C_dr = value(C_dr);
P_L = value(P_L);
size(P_en)
P_en = value(P_en);
P_loss = value(P_loss);
U;
