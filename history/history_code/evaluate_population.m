%% 计算目标函数值(目标函数包括包括年综合成本最小,电压偏差总量最小)
function [obj] = evaluate_population(LC_wt, LC_pv, t, P_wt, P_pv, P_en, P_loss, P_dr,U)  

% 初始化变量
N_s = 25;    
N_n = 33;

% 分时电价
d = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.45,0.56,0.56,0.45,0.45,0.4,0.56,0.4,0.45,0.45,0.45,0.4,0.4,0.2]; % 电价 单位元/1KWh



% 下层成本函数
% a.DG投资成本
% lamda_wt和lamda_pv为WT和PV的形状因子;alpha为贴现率;y是DG的使用寿命
% 现值系数:WT 0.06 PV 0.06 = 1/(1+alpha)^y;y=20 --贴现率:alpha ≈ 0.15(保留两位小数)
% lamda_wt = ((1 + alpha)^y * alpha) / ((1 + alpha)^y - 1);
% lamda_pv = ((1 + alpha)^y * alpha) / ((1 + alpha)^y - 1);
% lamda_wt = lamda_pv ≈ 0.16(保留两位小数)

% LC_wt(i)与LC_pv(i)为节点i安装的额定容量
% c_l_wt与c_l_pv为单位容量的投资成本
% C_l = lamda_wt * sum(c_l_wt .* LC_wt(i)) + lamda_pv * sum(c_l_pv .* LC_pv(i));
% C_l = 0.16 * sum(6200 .* LC_wt(i)) + 0.16 * sum(7000 .* LC_pv(i));
C_l = 24*(0.16 * sum(6200 .* LC_wt) + 0.16 * sum(7000 .* LC_pv)); % kw 和 kwh 单位换算

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
for s = 1:N_s
    for hour = 1:24
        C_en = C_en + d(1,hour) .* P_en(s,hour) .* t(s);
    end
end

% C_en = sum(sum(d .* P_en .* t, 2));


% d.网络损失费用
% 公式表示对所有场景的网络损失费用进行求和,其中包含损失系数(c_loss)、损失功率(P_loss)和场景s对应天数(t)
c_loss = 0.5;
C_loss = sum(c_loss .* P_loss .* t);

% e.实施需求响应的补偿成本
% P_dr(j,s)表示场景s下负载节点j处的可转移负载功率-可转移负载功率:配电网中能够通过需求响应(Demand Response, DR)机制在时间或空间上调整其用电行为的负荷功率
% c_dr表示每单位功率的可转移负载的补偿成本
C_dr = 0;
c_dr = 0.2;
for s = 1:N_s
    for j = 1:N_n
        C_dr = C_dr + c_dr * P_dr(s,j) * t(s);
    end
end

% 1.年综合成本
C_total = C_l + C_om + C_loss + C_en + C_dr;
F1 = C_total;


% 2.电压偏差总量
% abs(1 - U) 计算每个节点的电压偏差绝对值

F2 = 0;
for s = 1:N_s
    for i = 1:N_n
        for hour = 1:24
            F2 = F2 + abs(1 - U(s,i,hour)) * t(s);
        end
    end
end

% F2 = sum(sum(sum(abs(1 - V) .* t)));

obj = [F1,F2];
    
% 限制条件
    cv = [];
    % 1.DG装机容量约束条件 -安装容量最大值应作为已知参数出现-单位:kw
    for i = 1:N_n
        cv = [cv, 
          % P_wt_c(i) + P_pv_c(i) <= P_DG_c_max(i), 
          % P_wt_c(i) + P_pv_c(i) <= 800 %单位 kw
          % P_wt_c(i) + P_pv_c(i) >= 0
          % P_wt_c(i) <= P_wt_c_max(i),    
          % P_wt_c(i) >= 0,
          % P_pv_c(i) <= P_pv_c_max(i)
          % P_pv(i) >= 0

          % P_wt_c_max与P_pv_c_ma未给出具体值
          LC_wt(i) + LC_pv(i) <= 800
          LC_wt(i) + LC_pv(i) >= 0
          LC_wt(i) <= P_wt_c_max(i)
          LC_wt(i) >= 0
          LC_pv(i) <= P_pv_c_max(i)
          LC_pv(i) >= 0
         ];
    end
    cv = [cv; 
    sum(LC_wt) <= 1000
    sum(LC_pv) <= 700
    ]; 


    % 5.DG渗透率约束 45% 
    % xi 可再生能源在系统中的最大渗透率
    % P_z 是系统总负载-依靠其他参数来计算
    % 系统总有功和无功负荷为3715 + j2350 kVA-与IEEE33BW的负荷数据相冲突
    cv = [cv;
    % sum(LC_wt + LC_pv) <= xi*P_z
    sum(LC_wt + LC_pv) <= 0.45 * 3715
    ]; 


