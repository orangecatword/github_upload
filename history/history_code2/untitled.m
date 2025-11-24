% 加入动态重构部分

% 2025.10.3 问题是否应该考虑风光注入的无功-不需要考虑:要么是确定值依赖实例数据;要么是不确定值-放入鲁棒不确定性问题中;

% 2025.10.5
% 现有问题:发电机为什么是在耗电而不是发电;
% 为什么电流值计算公式之间存在差异？-电流值只在加入动态重构后的程序中出了问题
% 数据上的差异:接入有功负荷为IEEE33_VNOR的十分之一;同时由于计算规则的不同,无功负荷与IEEE33_VNOR中的值相差较大-通过分析原理,暂时不用专门为此修改


% 先不考虑:动态重构的时间;损失函数外的其他目标函数

% 1.简化一下程序:整个储能端的程序
% 1.5.修改一下生成图
% 2.改为动态重构 不同时刻开关不同

% 未来可能的扩展方向:
% 将每个储能的容量不一致
% 旧问题:是否要考虑风 光 储能的运行成本
clc
clear
tic
mpc = case33bw;

% 从节点流出的有功与无功功率
P = sdpvar(33, 24, 'full');
Q = sdpvar(33, 24, 'full');

% 支路功率
Pij = sdpvar(32, 24, 'full');
Qij = sdpvar(32, 24, 'full');
% Pij = sdpvar(37, 24, 'full');
% Qij = sdpvar(37, 24, 'full');

% 支路电流的平方
Iij = sdpvar(32, 24, 'full'); 
% Iij = sdpvar(37, 24, 'full');

% 支路电压的平方
U = sdpvar(33, 24, 'full'); 
% U0=1.0609; % 起始点电压-部分论文把首端电压调高 已起到抬高末端电压的作用

%% 定义动态重构的新增变量
% nb = 33;%节点数,根节点为33
% nl = 37;%支路数
% nc = 5;%联络开关数
% Zij=binvar(nl,1);%网架结构               
% Z0=[ones(nl-nc,1);zeros(nc,1)];%初始拓扑                      
% assign(Zij,Z0);      
% 
% % 发电机功率即为上级电网传送给节点1的功率
% Pg = sdpvar(nb,24);%发电机有功
% Qg = sdpvar(nb,24);%发电机无功
% Pgmax=[ones(1,24);zeros(32,24)];
% Qgmax=[ones(1,24);zeros(32,24)];
Umax=[1*1*ones(1,24);1.05*1.05*ones(32,24)];
Umin=[1*1*ones(1,24);0.95*0.95*ones(32,24)];
%%

branch = mpc.branch;
% 结论:直接从branch中获得的数据已经约等于标幺值了
r = branch(:,3);
x = branch(:,4);


% 获取负荷值-进一步复杂化,即节点负荷需求会随着时间而发生变化
pload = mpc.pload;% 负荷数据
pload_prim = mpc.pload_prim/(1000*10); % 10为基准值 最后得到的负荷为标幺值
qload_prim = mpc.qload_prim/(1000*10);
a = 3.715; % 单时段所有节点有功容量,MW
b = 2.3; % 单时段所有节点无功容量,MW
pload = pload/a;%得到各个时段与单时段容量的比例系数
qload = pload/b;%假设有功负荷曲线与无功负荷曲线相同
PL = pload_prim*pload;%得到33*24的负荷值,每一个时间段每个节点的负荷
QL = qload_prim*qload;

% 可转移负载功率 
P_dr_out = sdpvar(33, 24, 'full'); % 转移出的功率
P_dr_in = sdpvar(33, 24, 'full'); % 转移入的功率

% 可转移负荷状态
u_out = binvar(33, 24,'full');% 高峰转出
u_in = binvar(33, 24,'full');% 低峰转入

% 储能状态
Ies_c1 = binvar(1, 24,'full');% 充电状态
Ies_dc1 = binvar(1, 24,'full');% 放电状态
Ies_c2 = binvar(1, 24,'full');
Ies_dc2 = binvar(1, 24,'full');
Ies_c3 = binvar(1, 24,'full');
Ies_dc3 = binvar(1, 24,'full');

% 储能容量
Ees_max = 0.2; % 标幺值 实际为2MWh
% 储能设备剩余能量
Ees1 = sdpvar(1, 25,'full'); % 问题:加上24h完了后的下一时刻,根据时间-剩余容量图可知 25时刻很可能会跌出下限
Ees2 = sdpvar(1, 25,'full');
Ees3 = sdpvar(1, 25,'full');

% 储能功率
Pes_c1 = sdpvar(1, 24,'full');% 充电功率
Pes_dc1 = sdpvar(1, 24,'full');% 放电功率
Pes_c2 = sdpvar(1, 24,'full');
Pes_dc2 = sdpvar(1, 24,'full');
Pes_c3 = sdpvar(1, 24,'full');
Pes_dc3 = sdpvar(1, 24,'full');


% 支路链接情况
upstream = zeros(32,32);
dnstream = zeros(32,32);

% IEEE33BW版,符合现有的网络连接情况-同时根据后面的功率潮流约束公式 例如dnstream(1,18) = 1;即代表节点2与节点19有链接
for i = 1:32
    upstream(i,i) = 1;
end

for i=[1:16,18:20,22:23,25:31] % 即存在子节点的节点集
    dnstream(i,i+1) = 1;
end

% 分支支路
dnstream(1,18) = 1;
dnstream(2,22) = 1;
dnstream(5,25) = 1;
% upstream=zeros(33,32);%代表流入节点支路
% dnstream=zeros(33,32);%代表流出节点支路
% for i=1:32
%     upstream(i+1,i)=1;
% end
% 
% for i=[1:17,19:21,23:24,26:32]
%     dnstream(i,i)=1;
% end
% dnstream(2,18)=1;
% dnstream(3,22)=1;
% dnstream(6,25)=1;
%% 修改后的支路连接情况
% upstream=zeros(33,37);%代表流入节点支路
% dnstream=zeros(33,37);%代表流出节点支路
% for i=1:32
%     upstream(i+1,i)=1;
% end
% upstream(21,33)=1;%支路33为21-8支路，流入节点21
% upstream(15,34)=1;%支路34为15-9支路，流入节点15
% upstream(22,35)=1;%支路35为22-12支路，流入节点22
% upstream(33,36)=1;%支路36为33-18支路，流入节点33
% upstream(29,37)=1;%支路37为29-25支路，流入节点29
% 
% for i=[1:17,19:21,23:24,26:32]
%     dnstream(i,i)=1;
% end
% dnstream(2,18)=1;
% dnstream(3,22)=1;
% dnstream(6,25)=1;
% 
% %5条流入，对应5条流出
% dnstream(8,33)=1;
% dnstream(9,34)=1;
% dnstream(12,35)=1;
% dnstream(18,36)=1;
% dnstream(25,37)=1;
%%

% 在这里加上确定容量函数-放弃对风光容量的初始化-先放弃上层模型
% [LC_wt, LC_pv] = initialize_population(10);

% 确认风光的安装位置 风: 13 17 25 光: 4 7 27
Loc_pv_initial = [0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
Loc_wt_initial = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
% Ppv = (LC_pv(1,:)/1000)' * (mpc.pv/10); % 标幺值
% Pwt = (LC_wt(1,:)/1500)' * (mpc.wind/10);
Ppv = (Loc_pv_initial/3)' * (mpc.pv/10); % 标幺值
Pwt = (Loc_wt_initial/3)' * (mpc.wind/10);

% 确认储能的安装位置 (多安装几个储能试试) 6 16 29
Lc_pes1 = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Lc_pes2 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Lc_pes3 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];

% 约束条件
cv = [];
% 电压约束
cv = [cv, Umin <= U,U <= Umax];

%发电机功率约束           
% cv = [cv, -Pgmax <= Pg,Pg <= Pgmax,-Qgmax <= Qg,Qg <= Qgmax];
% %支路功率约束
% cv = [cv, -1*Zij*ones(1,24) <= Pij,Pij <= 1*Zij*ones(1,24);
%           -1*Zij*ones(1,24) <= Qij,Qij <= 1*Zij*ones(1,24)]; % 0.11这个数值是可变的
%支路电流约束
%cv = [cv, 0 <= Iij,Iij <= 1*Zij*ones(1,24)];
%% 储能约束
cv = [cv;
    Ies_c1 + Ies_dc1 <= 1; % 储能状态约束
    sum(Pes_c1) - sum(Pes_dc1) == 0; % 充放电的功率和最好一致
    Pes_c1 >= 0;
    Pes_c1 <= Ies_c1 * 0.02; % 储能容量的10% 
    Pes_dc1 >= 0;
    Pes_dc1 <= Ies_dc1 * 0.02;
    Ees1(1) == 0.5*0.2;
    Ees1 >= 0.2 * 0.2;
    Ees1 <= 0.8 * 0.2
    ]; 

cv = [cv;
    Ies_c2 + Ies_dc2 <= 1 % 储能状态约束
    sum(Pes_c2) - sum(Pes_dc2) == 0; % 充放电的功率和最好一致
    Pes_c2 >= 0;
    Pes_c2 <= Ies_c2 * 0.02; % 储能容量的10% 
    Pes_dc2 >= 0;
    Pes_dc2 <= Ies_dc2 * 0.02;
    Ees2(1) == 0.5*0.2;
    Ees2 >= 0.2 * 0.2;
    Ees2 <= 0.8 * 0.2
    ]; 

cv = [cv;
    Ies_c3 + Ies_dc3 <= 1 % 储能状态约束
    sum(Pes_c3) - sum(Pes_dc3) == 0; % 充放电的功率和最好一致
    Pes_c3 >= 0;
    Pes_c3 <= Ies_c3 * 0.02; % 储能容量的10% 
    Pes_dc3 >= 0;
    Pes_dc3 <= Ies_dc3 * 0.02;
    Ees3(1) == 0.5*0.2;
    Ees3 >= 0.2 * 0.2;
    Ees3 <= 0.8 * 0.2
    ]; 

cv = [cv;
    Ees1(2:25)==Ees1(1:24)+0.9*Pes_c1-1.1*Pes_dc1;
    Ees2(2:25)==Ees2(1:24)+0.9*Pes_c2-1.1*Pes_dc2;
    Ees3(2:25)==Ees3(1:24)+0.9*Pes_c3-1.1*Pes_dc3
    ];

%% 可转移负荷状态约束
cv = [cv;
    u_out + u_in <=1 % 可转移负荷状态约束
    sum(P_dr_out,2) - sum(P_dr_in,2) == 0;
    PL- P_dr_out + P_dr_in >= 0;
    P_dr_out >= 0;
    P_dr_out <= u_out*37.15*0.0001;
    P_dr_in >= 0;
    P_dr_in <= u_in*37.15*0.0001;
    P_dr_out(1,:) == 0;
    P_dr_in(1,:) == 0;
    ];


%% 二阶锥约束
% for i = 1:32
%     for t = 1:24
%         cv = [cv; cone([2*(Pij(i,t)); 2*(Qij(i,t));Iij(i,t) - U(branch(i,1),t)], Iij(i,t) + U(branch(i,1),t))]; 
%  % 我的想法:因为关于Iij的约束是二阶锥约束而非等式约束,所以Iij更像是在约束计算中的中间变量,最后的电流还得是功率的平方除以电压的平方
%     end
% end
%% 修改后的二阶锥约束
for i = 1:32
    for t = 1:24
        cv = [cv; cone([2*(Pij(i,t)); 2*(Qij(i,t));Iij(i,t) - U(branch(i,1),t)], Iij(i,t) + U(branch(i,1),t))]; 
    end
end
%%

%% 配电网辐射状结构约束
% cv = [cv, sum(Zij) == 32];                
% %提高运算速度            
% cv = [cv, Zij(1,1) == 1];         %这里可以确定的是，与主网相连的路一定是闭合的，所以可以提前置1，使运算更加快速 
% %分出12个section   将支路捆绑，然后设定其中多条串接的支路至多有一条断开      
% s1=3;s2=2;s3=1;s4=3;s5=3;s6=8;s7=4;s8=1;s9=2;s10=4;s11=4;s12=1;              
% cv = [cv, s1-sum(Zij(3:5)) <=1,s2-sum(Zij(6:7)) <=1];
% cv = [cv, s4-sum(Zij(9:11)) <=1,s5-sum(Zij(12:14)) <=1];
% cv = [cv, s6-sum(Zij(15:17))-sum(Zij(29:32))-sum(Zij(36))<=1];
% cv = [cv, s7-sum(Zij(18:20))-sum(Zij(2))<=1,s9-sum(Zij(21))-sum(Zij(35))<=1];
% cv = [cv, s10-sum(Zij(22:24))-sum(Zij(37))<=1,s11-sum(Zij(25:28))];

%%
cv = [cv;
    % 如果将P Q定义为从节点流出的有功和无功功率,那么熊壮壮公式更符合物理定义
    P(2:33,:) == upstream*Pij - upstream*(Iij .* (r*ones(1,24))) - dnstream * Pij; % 有无功功功率
    Q(2:33,:) == upstream*Qij - upstream*(Iij .* (x*ones(1,24))) - dnstream * Qij;
    % P == upstream*Pij - upstream*(Iij .* (r*ones(1,24))) - dnstream * Pij; 
    % Q == upstream*Qij - upstream*(Iij .* (x*ones(1,24))) - dnstream * Qij;
    U(branch(:,2),:) == U(branch(:,1),:) - 2.*r*ones(1,24).*Pij - 2.*x*ones(1,24).*Qij + (r.^2+x.^2)*ones(1,24).*Iij;  
    % 问题:在matpower的电压计算过程中 没有考虑r x前面的系数2-为什么？注意这个理论上的问题
    U(1,:) == 1;
    ];

%% 修改后的潮流约束
% 功率潮流约束
% cv = [cv;
% P + Pg == upstream*Pij - upstream*(Iij .* (r*ones(1,24))) - dnstream * Pij; 
% Q + Qg == upstream*Qij - upstream*(Iij .* (x*ones(1,24))) - dnstream * Qij;
% ];
% % 电压潮流约束
% m = 1.05*1.05 - 0.95*0.95;
% M = (ones(nl,24) - Zij*ones(1,24))*m;             
% cv = [cv, U(branch(:,1),:) - U(branch(:,2),:) <= M + 2*(r*ones(1,24)).*Pij + 2*(x*ones(1,24)).*Qij - ((r.^2 + x.^2))*ones(1,24).*Iij];
% cv = [cv, U(branch(:,1),:) - U(branch(:,2),:) >= -M + 2*(r*ones(1,24)).*Pij + 2*(x*ones(1,24)).*Qij - ((r.^2 + x.^2))*ones(1,24).*Iij];

%%
% 节点功率约束
cv = [cv;
    P == (PL - Ppv - Pwt) - P_dr_out + P_dr_in + Lc_pes1' * Pes_c1 - Lc_pes1' * Pes_dc1 +...  % 负荷需求转移走了 消纳的功率减小;反之负荷转入,PL增加-充电是在消耗功率
    + Lc_pes2' * Pes_c2 - Lc_pes2' * Pes_dc2 + Lc_pes3' * Pes_c3 - Lc_pes3' * Pes_dc3
    ];
cv = [cv;
    Q == QL
    ];


% 利用损失功率最小为优化目标 解出这些变量
I_custom = (Pij.^2 + Qij.^2)./U(branch(:,1),:); % 已经是标幺值,支路首端电流
P_loss = Iij.*(r*ones(1,24));
% 注:Iij与I_custom的全局相关性：r = 0.998（几乎完全线性相关）,比例关系稳定,所以P_loss的计算式也能反映实际的网损

% 目标函数-乘一年内发生该情况的天数
cost = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.45,0.56,0.56,0.45,0.45,0.4,0.56,0.4,0.45,0.45,0.45,0.4,0.4,0.2]; % 电价 单位元/1KWh-万元/10MWh(正好为标幺值单位)
% cost电价肯定是针对各个节点到用户
Puse = PL- P_dr_out + P_dr_in;
Puse_sum = sum(Puse,1); % 转换到1*24维
% 购电成本
f_cost = sum(Puse_sum.*cost);
% 损失成本
f_loss = 0.5 * sum(sum(P_loss)); % 猜测:0.5的单位也是元/kwh-万元/10MWh
% 可转移负荷成本(减少负荷的波动)
f_dr = 0.2 * sum(sum(P_dr_out)) + 0.8 * max(sum(P_dr_out,1));% 可转移负荷的惩罚 除了在总量上有惩罚,还应该在单时刻转移的量上惩罚
% 对照实验证明,需要在转移总量和转移峰值上均作出约束,才能实现预期效果 0.9为暂时确定的系数

f = f_cost + f_loss + f_dr ;

tic
ops = sdpsettings('solver', 'cplex', 'verbose', 0);
ops.cplex.preprocessing.presolve = 1; % 启用预处理
ops.cplex.workmem = 8192;  % 8GB
ops.cplex.timelimit = 60;  % 减少求解时间限制
ops.cplex.mip.tolerances.mipgap = 0.001;  % 放宽最优间隙
ops.cplex.parallel = 1;  % 启用并行计算
ops.cplex.mip.strategy.search = 1;  % 使用动态搜索
ops.cplex.mip.strategy.heuristicfreq = 100; % 调整启发式频率
optimize(cv, f, ops);
toc

U = value(U);
P = value(P);
Q = value(Q);
Pij = value(Pij);
Qij = value(Qij);
Iij = value(Iij);
I_custom = value(I_custom);

% Zij = value(Zij);
% Pg  = value(Pg);
% Qg  = value(Qg);

f = value(f);
P_loss = value(P_loss);

% 储能充放电功率
Pes_c1 = value(Pes_c1);
Pes_dc1 = value(Pes_dc1);
Ees1 = value(Ees1);
Pes_c2 = value(Pes_c2);
Pes_dc2 = value(Pes_dc2);
Ees2 = value(Ees2);
Pes_c3 = value(Pes_c3);
Pes_dc3 = value(Pes_dc3);
Ees3 = value(Ees3);

Puse = value(Puse); % 进行负荷转以后各个节点的负荷
P_dr_out = value(P_dr_out);
P_dr_in = value(P_dr_in);
% 把可转移负荷加入前后的时刻-功率图画出
PL_sum = sum(PL,1);
Puse_sum = sum(Puse,1);
% 总功率的时刻图 
P_sum = sum(P,1);

% 绘制PL_sum和Puse_sum折线图
figure;
hold on;
plot(1:24, PL_sum, 'b-o', 'LineWidth', 2, 'DisplayName', '原始负荷');
plot(1:24, Puse_sum, 'r-s', 'LineWidth', 2, 'DisplayName', '调整后负荷');
xlabel('时间 (小时)');
ylabel('总负荷 (标幺值)');
title('24小时负荷曲线对比');
legend('show');
grid on;
hold off;

figure
plot(1:24, U(33,:), 'g-o', 'LineWidth', 2, 'DisplayName', '节点33的电压');
% 原因:在16-20时刻 节点整体功率提高 导致电压在每个节点降的更多 
xlabel('时间 (小时)');
ylabel('电压 (标幺值)');
legend('show');
grid on;

figure
plot(1:25, Ees3, 'b-o', 'LineWidth', 2, 'DisplayName', '储能剩余能量');
xlabel('时间 (小时)');
ylabel('功率 (标幺值)');
legend('show');
grid on;

figure
plot(1:24, P_sum, 'b-o', 'LineWidth', 2, 'DisplayName', '功率');
xlabel('时间 (小时)');
ylabel('功率 (标幺值)');
legend('show');
grid on;
toc

