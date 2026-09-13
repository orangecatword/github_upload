function[objective] = function_up_BESS(x,idx_H,vals_H) 

%% 功能：将1×33向量分解为3个独热向量 + 非零值列表 下层
%% 1.设参
N = 10; % 样本数
objective = zeros(N,1);
f1 = zeros(5,1); % 定义子函数的目标函数
f2 = zeros(5,1);
f3 = zeros(5,1);
f4 = zeros(5,1);
f5 = zeros(5,1);

r_load1 = zeros(5,1); % 定义各个场景的故障恢复率
r_load2 = zeros(5,1);
r_load3 = zeros(5,1);
r_load4 = zeros(5,1);
r_load5 = zeros(5,1);

R_load1 = zeros(33,5); % 定义各个场景的节点故障恢复率
R_load2 = zeros(33,5);
R_load3 = zeros(33,5);
R_load4 = zeros(33,5);
R_load5 = zeros(33,5);

V_bias1 = zeros(33,5); % 定义各个场景节点的电压偏差
V_bias2 = zeros(33,5);
V_bias3 = zeros(33,5);
V_bias4 = zeros(33,5);
V_bias5 = zeros(33,5);

Rps1 = zeros(5,1); % 定义各个场景的供电可靠性
Rps2 = zeros(5,1);
Rps3 = zeros(5,1);
Rps4 = zeros(5,1);
Rps5 = zeros(5,1);

                          
% 建立20个样本
for i = 1:N

% x = [0	0	0	0	0	0.122372814151381	0	0.442246481586313	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.298942584145761	0.331494033628370	0	0	0	0	0	0.587744562955325	0	0	0];
% ===== 1️⃣ 找出非零元素位置与数值 =====
idx = find(x(i,:) ~= 0);      % 非零元素索引
vals = x(i,idx);              % 非零元素取值

% 2026.8.28得到的优化结果 objective r_load 
% idx = [8 23 30];
% vals = [0.5225     0.51139     0.46168];
% 2026.9.7 增加氢储能大容量以及缓启动性质后的优化结果
% idx = [6 8 24 30];
% vals = [0.11449 0.45127 0.48763 0.45135];
N_1 = [11; 23; 6; 3; 20];
N_2 = [29 2; 24 14; 22 7; 12 17; 18 26];
N_3 = [30 28 25; 23 14 3; 28 5 20; 30 10 19; 19 30 23];
N_4 = [21 3 22 9; 27 16 19 3; 25 27 22 26; 5 4 21 2; 8 30 16 9];
N_5 = [13 7 17 4 3; 15 29 26 12 21; 20 21 4 8 23; 24 26 28 27 18; 13 23 19 14 26];

tic
for j = 1:5
    [f1(j),r_load1(j),V_bias1(:,j),R_load1(:,j),Rps1(j)] = function_down_AOA(idx,vals,idx_H,vals_H,N_1(j,:));
    [f2(j),r_load2(j),V_bias2(:,j),R_load2(:,j),Rps2(j)] = function_down_AOA(idx,vals,idx_H,vals_H,N_2(j,:));
    [f3(j),r_load3(j),V_bias3(:,j),R_load3(:,j),Rps3(j)] = function_down_AOA(idx,vals,idx_H,vals_H,N_3(j,:));
    [f4(j),r_load4(j),V_bias4(:,j),R_load4(:,j),Rps4(j)] = function_down_AOA(idx,vals,idx_H,vals_H,N_4(j,:));
    [f5(j),r_load5(j),V_bias5(:,j),R_load5(:,j),Rps5(j)] = function_down_AOA(idx,vals,idx_H,vals_H,N_5(j,:));
end
toc

gamma = 0.08;   % 年现率
alpha = 20;     % 项目寿命，建议 20 年; 电池每10年更换1次
Ees_cost = 105; % 电池配置价格 单位:万美元/10MWh 
EH_cost   = 45; % 氢电池配置价格
% 初始投资成本
% C_inv = Ees_cost*sum(vals)+EH_cost*0.3;
C_inv = Ees_cost*sum(vals)+EH_cost*sum(vals_H); % 当前氢储能的容量为1
% 重置成本
C_rep = Ees_cost*sum(vals)/(1+gamma)^10+EH_cost*sum(vals_H)/(1+gamma)^10;
% 资本回收系数
CRF = (gamma*(1+gamma)^alpha)/((1+gamma)^alpha-1);
X = CRF*(C_inv+C_rep)/365; 

% 不同N-K情况概率占比 40.70% 33.20% 17.47% 6.67% 1.97%
% N-1: 11; 23; 6; 3; 20                                                             对应概率: 0.298	0.258	0.224	0.128	0.092
% N-2: 29-2; 24-14; 22-7; 12-17; 18-26                                              对应概率: 0.089	0.126	0.231	0.280	0.274
% N-3: 30-28-25; 23-14-3; 28-5-20; 30-10-19; 19-30-23                               对应概率: 0.227	0.154	0.289	0.200	0.130
% N-4: 21-3-22-9; 27-16-19-3; 25-27-22-26; 5-4-21-2; 8-30-16-9                      对应概率: 0.181	0.318	0.150	0.055	0.296
% N-5: 13-7-17-4-3; 15-29-26-12-21; 20-21-4-8-23; 24-26-28-27-18; 13-23-19-14-26    对应概率: 0.267	0.200	0.338	0.088	0.107
p1 = [0.298	0.258 0.224	0.128 0.092];
p2 = [0.089	0.126 0.231	0.280 0.274];
p3 = [0.227	0.154 0.289 0.200 0.130];
p4 = [0.181 0.318 0.150 0.055 0.296];
p5 = [0.267 0.200 0.338 0.088 0.107];

ob1 = p1*f1;
ob2 = p2*f2;
ob3 = p3*f3;
ob4 = p4*f4;  
ob5 = p5*f5;
objective(i,:) = 0.407*ob1 + 0.332*ob2 + 0.1747*ob3 + 0.0667*ob4 + 0.0197*ob5 + X;


r1 = p1*r_load1;
r2 = p2*r_load2;
r3 = p3*r_load3;
r4 = p4*r_load4;
r5 = p5*r_load5;
r_load = 0.407*r1 + 0.332*r2 + 0.1747*r3 + 0.0667*r4 + 0.0197*r5;

V1 = p1*V_bias1';
V2 = p2*V_bias2';
V3 = p3*V_bias3';
V4 = p4*V_bias4';
V5 = p5*V_bias5';
V_bias = 0.407*V1 + 0.332*V2 + 0.1747*V3 + 0.0667*V4 + 0.0197*V5;

R1 = p1*R_load1';
R2 = p2*R_load2';
R3 = p3*R_load3';
R4 = p4*R_load4';
R5 = p5*R_load5';
R_load = 0.407*R1 + 0.332*R2 + 0.1747*R3 + 0.0667*R4 + 0.0197*R5;

% 以往计算用户平均停电时间都是直接数数，现在不能继续采用此方法
% 计算供电可靠性Rps
P1 = p1*Rps1;
P2 = p2*Rps2;
P3 = p3*Rps3;
P4 = p4*Rps4;
P5 = p5*Rps5;
Rps = 0.407*P1 + 0.332*P2 + 0.1747*P3 + 0.0667*P4 + 0.0197*P5;

end




