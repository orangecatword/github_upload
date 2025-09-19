%% 1.根据IEEE33BW算出每个节点的P Q(牛拉法等)-牛拉法解潮流 找到可行解(最好找开源 官方数据)-从简单的程序开始加
% 可用matpower的case33bw直接计算P Q
% 通过matpower本身验证求解潮流计算公式的正确
% -去掉风光数据计算
    % 潮流约束-判断公式具体是否正确
    % 节点净负荷计算公式
    % 省去了拓扑结构,仅用两行公式总结
    for j=1:25
        P(:,:,j) = upstream*Pij(:,:,j) - upstream*(Iij(:,:,j).*(r*ones(1,24))) - dnstream*Pij(:,:,j);
    end
    % 考虑一下直接赋值 以及下面的公式修改
    Q = -Q_L;
    for j=1:25
        cv = [cv,Q(:,:,j) == upstream*Qij(:,:,j) - upstream*(Iij(:,:,j).*(x*ones(1,24))) - dnstream*Qij(:,:,j)];
    end
    % cv = [cv, Q == Q_L]; % 和大多数论文的约束公式对比,等式右侧整体加了个负号
    for j = 1:33
            if j ~= 33 % 若非平衡节点
                cv = [cv, P(j,:,:) == -P_L(j,:,:) + P_dr_in(j,:,:) - P_dr_out(j,:,:) + P_pv(j,:,:) + P_wt(j,:,:)]; % 和大多数论文的约束公式对比,等式右侧整体加了个负号
            else
                cv = [cv, P(j,:,:) == -P_L(j,:,:) + P_dr_in(j,:,:) - P_dr_out(j,:,:) + P_pv(j,:,:) + P_wt(j,:,:) + P_en(1,:,:)]; % 和大多数论文的约束公式对比,等式右侧整体加了个负号
            end
    end
    % 电压约束公式
    tic
    for j =1:25
        cv = [cv, U(branch(:,2),:,j) == U(branch(:,1),:,j) - 2*(r*ones(1,24)).*Pij(:,:,j) - 2*(x*ones(1,24)).*Qij(:,:,j) + ((r.^2+x.^2)*ones(1,24)).*Iij(:,:,j)];
        % cv = [cv, U(branch(:,1),:,j).*Iij(:,:,j) >= Pij(:,:,j).^2 + Qij(:,:,j).^2]; % 对熊壮壮论文(4-6)进行锥松弛
        % cv = [cv, 7 <= Pij(:,:,j).^2 + Qij(:,:,j).^2]; % 视在功率约束 7为标幺值
    end
        Pij_vec = reshape(Pij, [], 1);
        Qij_vec = reshape(Qij, [], 1);
        Iij_vec = reshape(Iij, [], 1);
        U_vec = reshape(U(branch(:,1),:,:), [], 1);
        %cv = [cv; norm([2*Pij_vec; 2*Qij_vec; Iij_vec-U_vec]) <= Iij_vec+U_vec];
        %cv = [cv; cone([Pij_vec; Qij_vec], 7)];
    toc

%% 2.P_en应该是被控量,P_dr是主动变量 加上Pg主动变量(带约束)

% 3.再去找优化算法,