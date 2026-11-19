function [theta1, theta2, dtheta1_max, dtheta2_max, phi_A, phi_B] = gamma_max_theta_speed(gamma1, gamma2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%GAMMA_MAX_THETA_SPEED 单位 gamma 速度圆下 theta 轴最大角速度分析（全解析）。
%
% 输入：
%   gamma1, gamma2 - 当前输出方向角（单位：度）
% 输出：
%   theta1, theta2 - 名义摆角（单位：度，取第一组逆解分支）
%   dtheta1_max    - theta1 的最大角速度（度/秒）
%   dtheta2_max    - theta2 的最大角速度（度/秒）
%   phi_A          - theta1 达到最大角速度时的 gamma 速度方向角（度）
%   phi_B          - theta2 达到最大角速度时的 gamma 速度方向角（度）
%
% 约束条件（写死）：gamma 速度矢量和为 1 度/秒，即
%   gamma1_dot^2 + gamma2_dot^2 = 1^2（度/秒）^2，
% gamma 速度方向扫一整圈：gamma_dot = (cos phi, sin phi) 度/秒，phi∈[0,360°)。
% 可选 model=coneik_model(alpha1,alpha2)，默认 20°/20°。
%
% 原理：速度映射 theta_dot = J^{-1}·gamma_dot 为线性映射，记 M = J^{-1}，则
%   theta1_dot = M(1,1)·cos phi + M(1,2)·sin phi
% 为 phi 的正弦型函数，其最大绝对值等于 M 第一行的模，
% 在 phi_A = atan2(M(1,2), M(1,1)) 方向达到；theta2_dot 同理（第二行）。
% 两个最大值独立计算、一般不在同一方向达到。几何意义：gamma 速度单位圆经 M
% 映射为 theta 速度椭圆，两个最大值即椭圆在 theta1、theta2 轴上的最大投影。
% 奇异位形（两圆锥相切，beta = aerfa1+aerfa2 或 |aerfa1-aerfa2|）J 退化，报错。

    if nargin < 3, model=[]; end
    gamma_speed = 1;   % gamma 速度矢量和（度/秒），写死

    [theta1, theta2, ~, ~] = inverse_motion_analysis(gamma1, gamma2, model);
    if cosd(rem(double(gamma1),360))==0
        error('gamma_max_theta_speed:ChartSingularity','Pitch/yaw velocity directions are undefined at an exact direction pole.');
    end
    J = jacobian_gamma_theta(theta1, theta2, model);
    if abs(det(J)) < 1e-12
        error('gamma_max_theta_speed: 当前位形接近奇异（两圆锥相切），速度无法唯一映射。');
    end
    M = gamma_speed * inv(J);

    dtheta1_max = hypot(M(1, 1), M(1, 2));
    phi_A = atan2d(M(1, 2), M(1, 1));
    dtheta2_max = hypot(M(2, 1), M(2, 2));
    phi_B = atan2d(M(2, 2), M(2, 1));
end
