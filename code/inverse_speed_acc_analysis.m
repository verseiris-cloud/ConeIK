function [theta1_dot_max, theta1_ddot_max, theta2_dot_max, theta2_ddot_max, info] = inverse_speed_acc_analysis(gamma1, gamma2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%INVERSE_SPEED_ACC_ANALYSIS Worst joint demand over both analytic branches.
%   Finite real scalar pitch/yaw inputs are in degrees; position-IK errors
%   propagate unchanged. For each branch, cover the entire task-velocity box
%   [-100,100]^2 deg/s and stationary task-acceleration box [-200,200]^2
%   deg/s^2. Outputs are the worst absolute joint requirements over the box
%   and both branches, in their original order. They are demands, not motor
%   capabilities. Acceleration assumes instantaneous ZERO joint velocity.
%
%   A rank-deficient branch cannot cover a two-dimensional task box. In that
%   case all four outputs are Inf, an infeasible-box sentinel, and optional
%   INFO.box_feasible is false. No projected, unmet target is used to compute
%   a finite requirement. Numerical breakdown at a full-rank branch raises a
%   distinct error. Requesting INFO does not change this numerical policy.

    if nargin < 3, model=[]; end
    [t11, t21, t12, t22, position_info] = inverse_motion_analysis(gamma1, gamma2, model);
    if cosd(rem(double(gamma1),360))==0
        error('inverse_speed_acc_analysis:ChartSingularity','Pitch/yaw motion capability is undefined at an exact direction pole.');
    end
    branches = [t11, t21; t12, t22];
    velocity_limits = [100; 100];
    acceleration_limits = [200; 200];
    corner_signs = [1, 1, -1, -1; 1, -1, 1, -1];
    velocity_corners = velocity_limits.*corner_signs;
    acceleration_corners = acceleration_limits.*corner_signs;
    targets = [velocity_corners, acceleration_corners];
    velocity_tolerance = 1e-8 + 1e-10*max(velocity_limits);
    acceleration_tolerance = 1e-7 + 1e-10*max(acceleration_limits);
    branch_velocity = zeros(2); branch_acceleration = zeros(2);
    box_feasible = true;
    if nargout > 4, branch_info = cell(1, 2); end

    for branch = 1:2
        [~, ~, J] = coneik_kinematics(branches(branch, :).',[0;0],model);
        state = coneik_jacobian_state(J);
        feasible = state.rank == 2;
        if feasible
            % A linear map's component extrema on a box occur at corners.
            % One solve reuses the same J for both velocity and acceleration.
            allocations = J \ targets;
            residual = J*allocations - targets;
            if any(~isfinite(allocations(:))) || any(~isfinite(residual(:))) || ...
                    any(max(abs(residual(:, 1:4)), [], 1) > velocity_tolerance) || ...
                    any(max(abs(residual(:, 5:8)), [], 1) > acceleration_tolerance)
                error('inverse_speed_acc_analysis:NumericalBreakdown', ...
                    'Full-rank branch %d failed finite-value or task-box residual checks.', branch);
            end
            branch_velocity(:, branch) = max(abs(allocations(:, 1:4)), [], 2);
            branch_acceleration(:, branch) = max(abs(allocations(:, 5:8)), [], 2);
            solver = 'direct_box_corners';
        else
            box_feasible = false;
            branch_velocity(:, branch) = Inf;
            branch_acceleration(:, branch) = Inf;
            % NaN explicitly marks allocations/residuals that were not solved;
            % it must not be interpreted as a numerical motion solution.
            allocations = nan(2, 8); residual = nan(2, 8);
            solver = 'infeasible_task_box';
        end
        if nargout > 4
            branch_info{branch} = struct('theta', branches(branch, :).', ...
                'jacobian', J, 'rank', state.rank, 'singular_values', state.s, ...
                'condition', state.condition, 'near_singular', state.near_singular, ...
                'rank_tolerance', state.rank_tolerance, 'box_feasible', feasible, ...
                'solver', solver, 'velocity_allocations', allocations(:, 1:4), ...
                'acceleration_allocations', allocations(:, 5:8), ...
                'velocity_residual', residual(:, 1:4), ...
                'acceleration_residual', residual(:, 5:8));
        end
    end

    if box_feasible
        velocity = max(branch_velocity, [], 2);
        acceleration = max(branch_acceleration, [], 2);
    else
        velocity = [Inf; Inf]; acceleration = [Inf; Inf];
    end
    theta1_dot_max = velocity(1); theta1_ddot_max = acceleration(1);
    theta2_dot_max = velocity(2); theta2_ddot_max = acceleration(2);
    if nargout > 4
        info = struct('box_feasible', box_feasible, 'position_info', position_info, ...
            'branches', [branch_info{:}], 'velocity_limits', velocity_limits, ...
            'acceleration_limits', acceleration_limits, ...
            'velocity_corners', velocity_corners, ...
            'acceleration_corners', acceleration_corners, ...
            'velocity_tolerance', velocity_tolerance, ...
            'acceleration_tolerance', acceleration_tolerance, ...
            'acceleration_assumption', 'instantaneous zero joint velocity', ...
            'current_velocity', [0; 0], 'damping', 0, ...
            'infeasible_sentinel', 'all four outputs are Inf');
    end
end
