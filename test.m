clc
clear
mpc = case33bw;
results = runpf('case33bw');
gen_power = results.gen(:,[2,3]); % 发电机的传输功率 即为支路1的有无功功率
bus_voltage = results.bus(:,8); % 第8列为电压幅值(pu)
bus_angle = results.bus(:,9);   % 第9列为电压相角(度)
actual_voltage = bus_voltage * mpc.bus(1,10); % 乘以基准电压(12.66kV)

% 实际值变标幺值还要再除以10(基准容量 mpc.baseMVA)
% 求解节点功率
% 节点注入功率(实际值)
bus_P = results.bus(:,3); % 有功功率(MW)
bus_Q = results.bus(:,4); % 无功功率(MVAr)

% 提取支路功率(实际值)
% 支路首端功率
branch_Pf = results.branch(:,14); % 有功功率(MW)
branch_Qf = results.branch(:,15); % 无功功率(MVAr)
% 支路末端功率
branch_Pt = results.branch(:,16); % 有功功率(MW)
branch_Qt = results.branch(:,17); % 无功功率(MVAr)

% 计算支路电流(标幺值)
branch = mpc.branch;
from_bus = branch(:,1);
to_bus = branch(:,2);

% 首端电流
If_pu = sqrt(branch_Pf.^2 + branch_Qf.^2) ./ bus_voltage(from_bus);
% 末端电流
It_pu = sqrt(branch_Pt.^2 + branch_Qt.^2) ./ bus_voltage(to_bus);

% 转换为实际值(A)
Ibase = mpc.baseMVA*1e6/(sqrt(3)*mpc.bus(1,10)*1e3); % 基准电流
If_actual = If_pu * Ibase;
It_actual = It_pu * Ibase;

% 输出结果
disp('支路电流计算结果:');
disp('支路编号 | 首端电流(A) | 末端电流(A)');
disp('----------------------------------');
for i = 1:length(If_actual)
    fprintf('%4d     | %9.2f   | %9.2f\n', i, If_actual(i), It_actual(i));
end
