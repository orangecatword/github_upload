%% 2025.10.2
%% 定义变量
% 定义支路功率方向变量
Pij_pos = sdpvar(37, 24, 'full'); % 正向功率
Pij_neg = sdpvar(37, 24, 'full'); % 反向功率
Qij_pos = sdpvar(37, 24, 'full'); % 正向功率
Qij_neg = sdpvar(37, 24, 'full'); % 反向功率
% 决定潮流方向的开关
pos = binvar(37, 24,'full');% 正向传输
neg = binvar(37, 24,'full');% 反向传输

% 实际支路功率为两者之差
Pij = Pij_pos - Pij_neg;
Qij = Qij_pos - Qij_neg;

switch_state = binvar(37, 24, 'full'); % 存在支路的开断情况
%% 支路连接情况
% 扩展upstream和dnstream矩阵到37×37
upstream = zeros(37,37);
dnstream = zeros(37,37);

% 保留原有32条支路的连接关系
for i = 1:32
    upstream(i,i) = 1;
end
for i=[1:16,18:20,22:23,25:31] % 即存在子节点的节点集
    dnstream(i,i+1) = 1;
end
% 分支支路
dnstream(1,18) = 1; % 节点2输出的是支路18
dnstream(2,22) = 1; % 节点3输出的是支路22
dnstream(5,25) = 1; % 节点6输出的是支路25

% 设置新增支路的连接关系（示例）
% 支路33
upstream(33,33) = 1; % 针对节点21 从33支路输入
dnstream(20,20) = 1; % 输出的是支路7和支路8
dnstream(20,21) = 1;
% 支路34
upstream(34,34) = 1; % 针对节点15 从34支路输入
dnstream(14,14) = 1; % 输出的是支路14和支路15
dnstream(14,15) = 1;
% 支路35
upstream(35,35) = 1; % 针对节点22 从35支路输入
dnstream(21,21) = 1; % 输出的是支路21
% 支路36
upstream(36,36) = 1; % 针对节点33 从36支路输入
dnstream(32,32) = 1; % 输出的是支路32
% 支路37
upstream(37,37) = 1; % 针对节点29 从37支路输入
dnstream(28,28) = 1; % 输出的是支路28和支路29
dnstream(28,29) = 1;
%% 约束条件
tic
% 潮流约束
for k = 1:37
    for t = 1:24  % 时间维度循环
        cv = [cv; 
            P(branchs(k,2),t) == upstream(k,:)*(switch_state(:,t).*(Pij_pos(:,t)-Pij_neg(:,t)))- ...
                  dnstream(k,:)*(switch_state(:,t).*(Pij_pos(:,t)-Pij_neg(:,t))) - ...
                  upstream(k,:)*(switch_state(:,t).*(Iij(:,t).*r))];
        cv = [cv;
            Q(branchs(:,2),t) == upstream(k,:)*(switch_state(:,t).*(Qij_pos(:,t)-Qij_neg(:,t))) - ...
                  dnstream(k,:)*(switch_state(:,t).*(Qij_pos(:,t)-Qij_neg(:,t))) - ...
                  upstream(k,:)*(switch_state(:,t).*(Iij(:,t).*x))];
        cv = [cv; implies(switch_state(k,t), ...
            U(branch(k,2),t) == U(branch(k,1),t) - 2*r(k)*(Pij_pos(k,t)-Pij_neg(k,t)) -...
            2*x(k)*ones(1,24).*(Qij_pos(k,t)-Qij_neg(k,t)) + (r(k).^2+x(k).^2).*Iij(k,t))
             ];
    end
end
toc

% 二阶锥约束
for i = 1:37
    for t = 1:24
        cv = [cv; cone([2*(Pij_pos(i,t) - Pij_neg(i,t)); 2*(Qij_pos(i,t) - Qij_neg(i,t));Iij(i,t) - U(branchs(i,1),t)], Iij(i,t) + U(branchs(i,1),t))]; 
 % 我的想法:因为关于Iij的约束是二阶锥约束而非等式约束,所以Iij更像是在约束计算中的中间变量,最后的电流还得是功率的平方除以电压的平方
    end
end


% 确保同一时刻只有一种方向有潮流功率
 cv = [cv;
     Pij_pos >= 0;
     Pij_pos <= pos; % 假设最大潮流为1
     Pij_neg >= 0;
     Pij_neg <= neg
     pos + neg <= 1;
     ];

% 对所有37条支路添加开关状态约束-大M法约束
M = 1; % 赋值给其足够大的常数，暂定为1
cv = [cv;
    Pij_pos <= M * switch_state;
    Pij_neg <= M * switch_state;
    ];

% 需要有电网辐射状结构的约束
% 确保37条支路中选择32条闭合（33节点-1个变电站）
cv = [cv; sum(switch_state,1) == 32];
% 确保所有潮流中所有节点的父节点只有一个
for t = 1:24
    for i = 1:33
        p = find(branch(:,2)==i);
        n = find(branch(:,1)==i);
    cv = [cv;sum(pos(p,t)) == 1;
             sum(neg(n,t)) == 1
         ];
    end
end

%% 新增的开关变换次数的目标函数
for t = 1:24
    switch_change(t) = sum(abs(switch_state(:,t+1)-switch_state(:,t)));
end
    f_switch_state = z*sum(switch_change); % 系数z表示对开关变化次数的严重性评估




    

%% 2025.10.4
% 程序新思路见文件CPLEX配电网重构单时段+多时段
%% 网络重构约束
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

%5条流入，对应5条流出
dnstream(8,33)=1;
dnstream(9,34)=1;
dnstream(12,35)=1;
dnstream(18,36)=1;
dnstream(25,37)=1;

Constraints = [Constraints, sum(Zij) == 32];                
%提高运算速度            
Constraints = [Constraints, Zij(1,1) == 1];         %这里可以确定的是，与主网相连的路一定是闭合的，所以可以提前置1，使运算更加快速 
%分出12个section   将支路捆绑，然后设定其中多条串接的支路至多有一条断开      
s1=3;s2=2;s3=1;s4=3;s5=3;s6=8;s7=4;s8=1;s9=2;s10=4;s11=4;s12=1;              
Constraints = [Constraints, s1-sum(Zij(3:5)) <=1,s2-sum(Zij(6:7)) <=1];
Constraints = [Constraints, s4-sum(Zij(9:11)) <=1,s5-sum(Zij(12:14)) <=1];
Constraints = [Constraints, s6-sum(Zij(15:17))-sum(Zij(29:32))-sum(Zij(36))<=1];
Constraints = [Constraints, s7-sum(Zij(18:20))-sum(Zij(2))<=1,s9-sum(Zij(21))-sum(Zij(35))<=1];
Constraints = [Constraints, s10-sum(Zij(22:24))-sum(Zij(37))<=1,s11-sum(Zij(25:28))];

%欧姆定律约束
m = 1.05*1.05 - 0.95*0.95;
M = (ones(nl,T) - Zij*ones(1,T))*m;             
Constraints = [Constraints, V(branch(:,1),:) - V(branch(:,2),:) <= M + 2*(r*ones(1,T)).*Pij + 2*(x*ones(1,T)).*Qij - ((r.^2 + x.^2))*ones(1,T).*Iij];
Constraints = [Constraints, V(branch(:,1),:) - V(branch(:,2),:) >= -M + 2*(r*ones(1,T)).*Pij + 2*(x*ones(1,T)).*Qij - ((r.^2 + x.^2))*ones(1,T).*Iij];


%% 2025.10.4
% 预计算常用表达式
% upstream_switch = upstream .* switch_state';
% dnstream_switch = dnstream .* switch_state';
% Pij_diff = Pij_pos - Pij_neg;
% Qij_diff = Qij_pos - Qij_neg;
% 
% % 向量化计算
% for t = 1:24
%     cv = [cv;
%         P(branchs(:,2),t) == upstream_switch * Pij_diff(:,t) - dnstream_switch * Pij_diff(:,t) - upstream_switch * (Iij(:,t).*r);
%         Q(branchs(:,2),t) == upstream_switch * Qij_diff(:,t) - dnstream_switch * Qij_diff(:,t) - upstream_switch * (Iij(:,t).*x)];
% 
%     % 保持implies逻辑
    for k = 1:37
        cv = [cv; implies(switch_state(k,t), ...
            U(branch(k,2),t) == U(branch(k,1),t) - 2*r(k)*Pij_diff(k,t) - ...
            2*x(k)*Qij_diff(k,t) + (r(k)^2+x(k)^2)*Iij(k,t))];
    end
% end

%% 2025.10.3
% for t = 1:24
%     active_P = switch_state(:,t).*(Pij_pos(:,t)-Pij_neg(:,t));
%     active_Q = switch_state(:,t).*(Qij_pos(:,t)-Qij_neg(:,t));
%     active_I = switch_state(:,t).*Iij(:,t);
% 
%     for k = 1:37
%         cv = [cv; 
%             P(branchs(k,2),t) == upstream(k,:)*active_P - dnstream(k,:)*active_P - upstream(k,:)*(active_I.*r);
%             Q(branchs(k,2),t) == upstream(k,:)*active_Q - dnstream(k,:)*active_Q - upstream(k,:)*(active_I.*x);
%         ];
% 
%         if switch_state(k,t)
%             cv = [cv; 
%                 U(branch(k,2),t) == U(branch(k,1),t) - 2*r(k)*(Pij_pos(k,t)-Pij_neg(k,t)) -...
%                 2*x(k)*(Qij_pos(k,t)-Qij_neg(k,t)) + (r(k)^2+x(k)^2)*Iij(k,t)
%             ];
%         end
%     end
% end

%% 2025.10.1
% cv = [cv; implies(switch_state(k,:), 
%        U(branch(k,2),:) == U(branch(k,1),:) - 2*r(k)*(Pij_pos(k,:)-Pij_neg(k,:)) - ...)];
% cv = [cv; P(branchs(i,2),:) == sum(upstream(:,i).*switch_state(:,t)'.*(Pij_pos-Pij_neg),1) - ...];

