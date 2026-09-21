%% 基于论文-Optimal planning of hybrid hydrogen and battery energy storage for resilience enhancement using bi-layer decomposition algorithm 

% ---------- 决策变量 ----------
% N_H：现场氢站数量；T：调度时域长度；delta_t：单个时段长度(h)

Pel = sdpvar(N_H, T, 'full');       % 电解槽耗电功率 P^E (MW)
Pfc = sdpvar(N_H, T, 'full');       % 燃料电池发电功率 P^F (MW)
EH  = sdpvar(N_H, T+1, 'full');     % 储氢量 E^H (MWh_H2)

u_el = binvar(N_H, T, 'full');      % 电解槽运行状态
u_fc = binvar(N_H, T, 'full');      % 燃料电池运行状态

EHS = sdpvar(N_H, T, 'full');   % 售出的氢能量，单位应与 EH 一致


% ---------- 参数 ----------
% Pel_max : N_H×1，电解槽额定耗电功率
% Pfc_max : N_H×1，燃料电池额定发电功率
% EH_max  : N_H×1，储氢罐额定容量
% eta_el  : N_H×1，电解制氢效率
% eta_fc  : N_H×1，燃料电池发电效率
% rho0    : 初始储氢比例，例如 0.5
% rho_min : 最低储氢比例，例如 0 或 0.2
% rho_max : 最高储氢比例，通常为 1
% sigma_H : 每小时储氢损耗率；若忽略损耗，取 0

eta_el_mat = eta_el * ones(1, T);
eta_fc_mat = eta_fc * ones(1, T);
sigma_H_mat = (1 - sigma_H * delta_t) * ones(N_H, T);

% ---------- 现场氢储能约束 ----------
Constraints = [Constraints;

    % 电解槽与燃料电池出力范围
    0 <= Pel <= Pel_max * ones(1, T);
    0 <= Pfc <= Pfc_max * ones(1, T);

    % 启停状态与功率绑定
    Pel <= (Pel_max * ones(1, T)) .* u_el;
    Pfc <= (Pfc_max * ones(1, T)) .* u_fc;

    % 同一现场氢站同一时段不能同时制氢、发电
    u_el + u_fc <= 1;

    % 初始储氢量
    EH(:, 1) == rho0 * EH_max;

    % 储氢容量边界
    rho_min * EH_max * ones(1, T+1) <= EH <= ...
    rho_max * EH_max * ones(1, T+1);

    % 储氢动态：制氢增加氢能，燃料电池发电消耗氢能
    EH(:, 2:T+1) == sigma_H_mat .* EH(:, 1:T) ...
                  + eta_el_mat .* Pel * delta_t ...
                  - (Pfc ./ eta_fc_mat) * delta_t...
                  - EHS;
];
