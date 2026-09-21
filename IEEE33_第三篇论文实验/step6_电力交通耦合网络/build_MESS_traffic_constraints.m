function [Constraints, y_mess, x_mess, travel_step] = ...
    build_MESS_traffic_constraints(traffic,T,Delta_t,Nmess,init_node,road_available,t0_scale)

% BUILD_MESS_TRAFFIC_CONSTRAINTS
% 建立MESS在Sioux Falls交通网络中的移动约束
%
% 输入：
% traffic        Sioux Falls交通网络结构体
% T              总时间步数
% Delta_t        每个时间步对应的实际时间
% Nmess          MESS数量
% init_node      每台MESS初始交通节点
% road_available 道路可用状态，大小为Nlink×T
% t0_scale       自由流时间缩放系数
%
% 输出：
% Constraints    YALMIP约束集合
% y_mess{m}(i,t) 第m台MESS在时刻t是否位于节点i
% x_mess{m}(l,t) 第m台MESS在时刻t是否进入道路l
% travel_step(l) 第l条道路所需的离散时间步数

%% 1. 基础参数
branch_t = traffic.branch;
Nnode = traffic.n_node;
Nlink = traffic.n_link;

from = branch_t(:,traffic.col.from);
to = branch_t(:,traffic.col.to);
t0 = branch_t(:,traffic.col.t0);

%% 2. 输入检查
if length(init_node) ~= Nmess
    error('init_node长度必须等于Nmess');
end

if size(road_available,1) ~= Nlink
    error('road_available行数必须等于交通道路数');
end

if size(road_available,2) == 1
    road_available = repmat(road_available,1,T);
end

if size(road_available,2) ~= T
    error('road_available列数必须等于T');
end

%% 3. 道路实际行驶时间
travel_time = t0_scale .* t0;

%% 4. 转换为离散时间步
travel_step = ceil(travel_time ./ Delta_t);
travel_step = max(travel_step,ones(Nlink,1));

%% 5. 定义YALMIP变量
y_mess = cell(Nmess,1);
x_mess = cell(Nmess,1);

for m = 1:Nmess
    y_mess{m} = binvar(Nnode,T,'full');
    x_mess{m} = binvar(Nlink,T,'full');
end

%% 6. 初始化约束
Constraints = [];

%% 7. MESS初始位置约束
for m = 1:Nmess
    y0 = zeros(Nnode,1);
    y0(init_node(m)) = 1;
    Constraints = [Constraints, y_mess{m}(:,1) == y0];
end

%% 8. 每台MESS同一时刻最多位于一个节点
% 当MESS位于道路中时，sum(y_mess(:,t))=0
for m = 1:Nmess
    for t = 1:T
        Constraints = [Constraints, sum(y_mess{m}(:,t)) <= 1];
    end
end

%% 9. 道路故障约束
% road_available(l,t)=0时，禁止进入道路l
for m = 1:Nmess
    Constraints = [Constraints, x_mess{m}(:) <= road_available(:)];
end

%% 10. 只有位于某节点时，MESS才能从该节点出发
% sum_{l∈out(i)} x(m,l,t) <= y(m,i,t)
for m = 1:Nmess
    for t = 1:T
        for i = 1:Nnode
            out_idx = find(from == i);
            if ~isempty(out_idx)
                Constraints = [Constraints, ...
                    sum(x_mess{m}(out_idx,t)) <= y_mess{m}(i,t)];
            end
        end
    end
end

%% 11. 禁止在仿真时间结束后才到达
% 要求：t + travel_step(l) <= T
for m = 1:Nmess
    for l = 1:Nlink
        latest_depart = T - travel_step(l);
        if latest_depart < T
            first_forbidden = max(1,latest_depart + 1);
            Constraints = [Constraints, ...
                x_mess{m}(l,first_forbidden:T) == 0];
        end
    end
end

%% 12. MESS节点状态动态约束
% y(i,t+1)
% = y(i,t)
% - 当前从节点i离开的MESS
% + 当前到达节点i的MESS

for m = 1:Nmess
    for t = 1:T-1
        for i = 1:Nnode

            % 从节点i离开的MESS
            out_idx = find(from == i);
            outflow = 0;
            if ~isempty(out_idx)
                outflow = sum(x_mess{m}(out_idx,t));
            end

            % 在t+1时刻到达节点i的MESS
            in_idx = find(to == i);
            inflow = 0;

            for kk = 1:length(in_idx)
                l = in_idx(kk);
                departure_time = t + 1 - travel_step(l);

                if departure_time >= 1
                    inflow = inflow + x_mess{m}(l,departure_time);
                end
            end

            % 节点状态更新
            Constraints = [Constraints, ...
                y_mess{m}(i,t+1) == ...
                y_mess{m}(i,t) - outflow + inflow];
        end
    end
end

end