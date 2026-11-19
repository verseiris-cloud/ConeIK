function [PA, Opt, s1, s2] = circle_from_vector(n, alpha_deg, theta_deg)
% CIRCLE_FROM_VECTOR 根据输入向量、夹角参数和摆角计算圆锥圆的几何量。
%
% 输入：
%   n          - 三维非零向量（行向量或列向量均可）
%   alpha_deg  - 夹角参数（单位：度），即圆锥半顶角
%   theta_deg  - 摆角（单位：度），决定圆上点的位置
%
% 输出：
%   PA  - 圆心三维坐标（列向量），对应文档中的 P_A
%   Opt - 圆上目标点三维坐标（列向量），对应文档中的 Opt
%   s1  - 圆所在平面内第一个正交基向量，长度为 R_A
%   s2  - 圆所在平面内第二个正交基向量，长度为 R_A，且 s2 与 s1 垂直
%
% 计算步骤：
%   1. 将输入向量 n 单位化得到 u。
%   2. 以 u 为轴、半顶角 alpha_deg 的圆锥与单位球面的交线是一个圆，
%      圆心 PA = cosd(alpha_deg) * u，半径 RA = sind(alpha_deg)。
%   3. 选取特殊方向 n0 = [1,0,0]^T，对应正交基为 s1_0 = [0,RA,0]^T、
%      s2_0 = [0,0,RA]^T（均垂直于 n0 且长度为 RA）。
%      通过 Rodrigues 旋转，将 n0 旋转到 u，同步旋转 s1_0、s2_0，
%      得到当前输入向量对应的正交基 s1、s2。
%   4. 圆上目标点 Opt = PA + cosd(theta_deg)*s1 + sind(theta_deg)*s2。

    n = n(:);                         % 统一转为列向量
    nv = norm(n);
    if nv < eps
        error('circle_from_vector: 输入向量不能为零向量。');
    end
    u = n / nv;

    % 角度制转弧度制
    alpha = deg2rad(alpha_deg);
    theta = deg2rad(theta_deg);

    ca = cos(alpha);
    sa = sin(alpha);
    PA  = ca * u;                     % 圆心坐标
    RA  = sa;                         % 圆半径

    % 构造旋转：将标准方向 n0 = [1;0;0] 旋转到 u
    n0 = [1; 0; 0];
    k = cross(n0, u);
    nk = norm(k);
    if nk < eps
        % u 与 n0 平行或反平行
        if dot(n0, u) > 0
            s1 = [0; RA; 0];
            s2 = [0; 0; RA];
        else
            % 反平行时相当于绕 y 轴旋转 180 度
            s1 = [0;  RA;  0];
            s2 = [0;  0; -RA];
        end
    else
        k = k / nk;
        angle = atan2(nk, dot(n0, u)); % 从 n0 到 u 的旋转角
        R = rodrigues_matrix(k, angle);
        s1 = R * [0; RA; 0];
        s2 = R * [0; 0; RA];
    end

    % 计算圆上目标点
    Opt = PA + cos(theta) * s1 + sin(theta) * s2;
end

function R = rodrigues_matrix(k, angle)
%RODRIGUES_MATRIX 罗德里格斯旋转矩阵
%   k 为单位旋转轴（列向量），angle 为旋转角（弧度）。
    K = [0    -k(3)  k(2);
         k(3)  0    -k(1);
        -k(2)  k(1)  0   ];
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K * K);
end
