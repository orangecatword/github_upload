% 接下来参数对接与归一化-实现了价格部分的调整
% 问题2:如何计算某节点i在某场景s下的光伏/风电功率 已知不同时刻节点在某场景下的光伏/风电功率 -可以通过求取24h的平均值来解决问题
% 问题3:单个节点的光伏与风电容量限制分别为? -论文未给出, 后续也许可以自己设置
% 问题4:某个场景下节点的实时电压如何计算 - 不同场景的不同节点求最值即可


% 参考文献:英_配电网中多目标分布式发电分层最优规划:改进的白鲸优化算法

% 论文在案例分析中给出的参数:系统总有功和无功负荷为3715 + j2350 kVA;
% IEEE 33在IEEE3BW.m中,第一行是母线(即节点1)与节点2的线路数据    -在论文《英_采用多目标麻雀搜索算法建立主动配电网动态重构集成优化》中给出了IEEE33线路的视在功率
% 系统参考电压为12.66 kV,作为电压归一化分母
% WT的安装节点位置至多为13、17、25;不同的优化算法可能会选择其中的两个点;系统安装的WT最大总容量为1000kw 
% PV的安装节点数为4、7、27;不同的优化算法可能会选择其中的1-2个点;每个节点的DG装机容量上限为800kw, DG总容量上限为1700kw

% 投资成本(元/kwh):WT 6200 PV 7000
% 运营和维护成本(元/kwh):WT 0.2 PV 0.2
% 最大使用寿命(年):WT 20 PV 20
% 现值系数:WT 0.06 PV 0.06

% 分时电价-已放入程序

% 系统电压限制(p.u): 0.95~1.05
% 分支电路容量限制(MVA) 0~7
% 单位需求响应补偿成本(元)0.2元
% 节点最大可转移负荷(kW)37.15KW
% 单位网络损耗成本(元)0.5元
% 最大可再生能源渗透率 45%
% 其中，除了分布式电源外,暂无其他的调控设备

clc
clear
warning off
tic
global km kn
kn=4;%循环次数
%% 初始化参数
mpc = IEEE33BW; %IEEE33标准节点系统

% wind和pv都是不确定性变量d的初始值
wind = mpc.wind; % 风电数据
wind = wind/3; % 含3个风机
wind = ones(3,1).*wind;

pv = mpc.pv;%光伏数据,含3个
pv = pv/3;
pv = ones(3,1).*pv;

% 负荷定义不确定性变量d的初始值
pload = mpc.pload;% 负荷数据
pload_prim = mpc.pload_prim/1000;%化为标幺值
qload_prim = mpc.qload_prim/1000;
a = 3.715;%单时段所有节点有功容量,MW
b = 2.3;%单时段所有节点无功容量,MW
pload = pload/a;%得到各个时段与单时段容量的比例系数
qload = pload/b;%假设有功负荷曲线与无功负荷曲线相同
pload = pload_prim*pload;%得到33*24的负荷值,每一个时间段每个节点的负荷
qload = qload_prim*qload;
pload = pload(1:32,:);   
qload = qload(1:32,:);



branch = mpc.branch;
branch(:,3) = branch(:,3)*1/(12.66^2);%求阻抗标幺值 论文中电压压基准值选择为12.66kV
r = real(branch(:,3));
x = imag(branch(:,3));

% 求导纳的实部与虚部
g = r./(r.^2 + x.^2);
b = -x./(r.^2 + x.^2);

T = 24; %时段数为24小时
nb = 33;%节点数
nl = 32;%支路数

nwt = 3; % 3个风机
npv = 3; % 3个光伏

% 先在main中初始化变量-最终都是要将这些内容放到上下层模型中
% LC 代指位置和容量
LC_wt = sdpvar(33, 1); % 风机对应位置的容量
LC_pv = sdpvar(33, 1); % 光伏对应位置的容量

% 节点的风电与光伏额定装机容量
% P_wt_c = sdpvar(33, 1); 即为LC_wt
% P_pv_c = sdpvar(33, 1); 即为LC_pv
% P_DG_c = sdpvar(33, 1); LC_wt + LC_pv

% 允许在节点安装的DG 风力 光伏最大容量
P_DG_c_max = sdpvar(33, 1);
P_wt_c_max = sdpvar(33, 1);
P_pv_c_max = sdpvar(33, 1);

% 节点的风电与光伏出力(某个场景的某个节点)
P_wt = sdpvar(25, 33); 
P_pv = sdpvar(25, 33);
P_DG = sdpvar(25, 33);

% 分布式电源的有功出力上限(某个场景的某个节点)
P_DG_max = sdpvar(25, 33);

% 实时电压
V = sdpvar(25, 33, 24);

% 
P = sdpvar(25, 33); % 场景s时流入节点i的有功功率和无功功率
Q = sdpvar(25, 33);
U = sdpvar(25, 33); % 场景s时流入节点i的电压

% 支路视在功率
Sij = sdpvar(32, 1); 

% 分布式电源的切削率
omega_DG = sdpvar(33, 1);
omega_DG_max = sdpvar(33, 1);

% 可转移负载功率
P_dr = sdpvar(25, 33);

% 系统的实际价格
d = sdpvar(25, 24);

% 上级电网购买电力的有功功率
P_en = sdpvar(25, 24);

% 损失功率
P_loss = sdpvar(1, 25);

% 不同场景的运行天数
t = sdpvar(1, 25);

% 不同场景出现的概率
pk = sdpvar(1, 25);
% C_res代表聚类中心-大小设置无十足把握
% 聚类中心代表其中一个节点(光伏/风电)24小时的出力情况
C_res = sdpvar(25, 48); % 场景总数为25,风光各占24h

% 获取光伏风电的不确定性在典型场景下的出力情况
% pk:不同场景出现的概率
% 不同场景下 光伏和风力的出力值
[C_res,pk] = datap; % 场景总数为25

% t =  365 * pk 如果要是天数,还需考虑四舍五入- pk 为概率,t为天数 能不能通过加一个约束来解决 天数和为365 同时 t ≈ 365 * pk(不能差超过1天)

% 导入的的因素是场景下光伏 风力出力值
[new_population, new_obj] = up_configuration(C_res,pk,g,b);


% 结果
V = value(V);%电压的平方
I = value(I);%电流的平方
P = value(P);%线路有功
Q = value(Q);%线路无功
p_wt = value(p_wt);%风机有功
p_pv = value(p_pv);%光伏有功
Pin = value(Pin);
Qin = value(Qin);
P_wt = value(P_wt);
P_pv = value(P_pv);



%画图



% 图2:三维电压分布图-展示各节点在24小时内的电压幅值时空分布
figure;
yy=1:24;
xn=1:33;
% mesh:生成三维网格曲面图,横轴为时间(1-24小时),纵轴为节点编号(1-33),竖轴为电压标幺值(V)
mesh(yy,xn,V(:,:,end)),xlabel('时间'),ylabel('节点'),zlabel('电压(p.u.)') % V(:,:,end):取最后一次迭代的电压数据。
title('电压');

% 图3:线路有功功率分布图
figure;
yy=1:24;
xn=1:32;
mesh(yy,xn,P(:,:,end)),xlabel('时间'),ylabel('节点'),zlabel('线路有功')

% 图4:功率平衡堆叠图-可视化24小时内各电源出力与负荷的平衡关系
figure;
% 堆叠柱状图:bar(yyf, 'stack')绘制储能充电功率(负值),yyz包含储能放电、风电(3台)、光伏(2台)和购电功率
yyf=[-p_ch(1,:,end);-p_ch(2,:,end)]';
bar(yyf,'stack');
yyz=[p_dis(1,:,end);p_dis(2,:,end);p_wt(1,:,end);p_wt(2,:,end);p_wt(3,:,end);p_pv(1,:,end);p_pv(2,:,end);Pg(end,:,end)]';
hold on
bar(yyz,'stack');
plot(sum(pload),'b--*','LineWidth',1.5)
legend('储能1充电','储能2充电','储能1放电','储能2放电','风电1','风电2','风电3','光伏1','光伏2','购电','有功负荷');
grid on
xlabel('时间');
ylabel('功率');

