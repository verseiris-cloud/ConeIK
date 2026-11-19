function [gamma1_dot, gamma1_ddot, gamma2_dot, gamma2_ddot, info] = speed_acc_analysis(theta1, theta2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%SPEED_ACC_ANALYSIS Independent task-axis capability with joint box limits.
%   Inputs are finite real scalar NOMINAL angles, in degrees. The first four
%   outputs retain their order and are nonnegative speeds (deg/s) and
%   accelerations (deg/s^2). Each task axis moves alone while the other task
%   component stays zero. Joint limits are [300;700] deg/s and
%   [1500;3600] deg/s^2.
%
%   Acceleration capability is instantaneous at ZERO joint velocity. At a
%   nonzero velocity, the curvature term changes the feasible acceleration
%   set, so these values must not be used as moving-state acceleration limits.
%   A full-rank J gives the unique allocation for each task direction. At
%   rank one, a supported direction uses the exact joint-box support value,
%   including nullspace allocation; an unsupported direction has capability 0.
%   Optional INFO reports rank, conditioning, allocations and residuals.

    if nargin < 3, model=[]; end
    if ~(isnumeric(theta1) && isreal(theta1) && isscalar(theta1) && ...
            isfinite(theta1) && isnumeric(theta2) && isreal(theta2) && ...
            isscalar(theta2) && isfinite(theta2))
        error('speed_acc_analysis:InvalidInput', ...
            'theta1 and theta2 must be finite real numeric scalars in degrees.');
    end
    velocity_limits = [300; 700];
    acceleration_limits = [1500; 3600];
    [~, ~, J] = coneik_kinematics([double(theta1); double(theta2)],[0;0],model);
    state = coneik_jacobian_state(J);
    targets = eye(2);
    supported = false(1, 2);
    velocity = zeros(1, 2); acceleration = zeros(1, 2);
    velocity_allocation = zeros(2); acceleration_allocation = zeros(2);

    if state.rank == 2
        unit_allocation = J \ targets;
        if any(~isfinite(unit_allocation(:)))
            error('speed_acc_analysis:NumericalBreakdown', ...
                'The full-rank axis allocation produced a nonfinite result.');
        end
        for axis = 1:2
            magnitudes = abs(unit_allocation(:, axis));
            % An exact zero coefficient imposes no restriction. Do not
            % threshold small coefficients or alter near-singular geometry.
            velocity(axis) = min(velocity_limits ./ magnitudes);
            acceleration(axis) = min(acceleration_limits ./ magnitudes);
            velocity_allocation(:, axis) = unit_allocation(:, axis)*velocity(axis);
            acceleration_allocation(:, axis) = unit_allocation(:, axis)*acceleration(axis);
        end
        supported(:) = true;
        solver = 'direct';
    elseif state.rank == 1
        range_basis = state.U(:, state.active);
        % This is only a floating-point range-membership tolerance; the
        % numerical rank threshold is supplied by coneik_jacobian_state.
        direction_tolerance = 64*eps;
        for axis = 1:2
            target = targets(:, axis);
            supported(axis) = norm(target-range_basis*(range_basis.'*target), inf) ...
                <= direction_tolerance;
            if supported(axis)
                coefficients = target.'*J;
                velocity(axis) = sum(abs(coefficients).*velocity_limits.');
                acceleration(axis) = sum(abs(coefficients).*acceleration_limits.');
                velocity_allocation(:, axis) = sign(coefficients(:)).*velocity_limits;
                acceleration_allocation(:, axis) = sign(coefficients(:)).*acceleration_limits;
            end
        end
        solver = 'rank_one_box_support';
    else
        solver = 'zero_rank';
    end

    velocity_residual = J*velocity_allocation - diag(velocity);
    acceleration_residual = J*acceleration_allocation - diag(acceleration);
    velocity_tolerance = 1e-8 + 1e-10*velocity;
    acceleration_tolerance = 1e-7 + 1e-10*acceleration;
    values = [velocity, acceleration, velocity_allocation(:).', ...
        acceleration_allocation(:).', velocity_residual(:).', acceleration_residual(:).'];
    if any(~isfinite(values)) || ...
            any(max(abs(velocity_residual), [], 1) > velocity_tolerance) || ...
            any(max(abs(acceleration_residual), [], 1) > acceleration_tolerance)
        error('speed_acc_analysis:NumericalBreakdown', ...
            'The computed capability failed finite-value or achieved-motion checks.');
    end

    gamma1_dot = velocity(1); gamma1_ddot = acceleration(1);
    gamma2_dot = velocity(2); gamma2_ddot = acceleration(2);
    if nargout > 4
        info = struct('jacobian', J, 'rank', state.rank, ...
            'singular_values', state.s, 'condition', state.condition, ...
            'near_singular', state.near_singular, ...
            'rank_tolerance', state.rank_tolerance, 'solver', solver, ...
            'supported_directions', supported, 'velocity_limits', velocity_limits, ...
            'acceleration_limits', acceleration_limits, ...
            'velocity_allocation', velocity_allocation, ...
            'acceleration_allocation', acceleration_allocation, ...
            'velocity_residual', velocity_residual, ...
            'acceleration_residual', acceleration_residual, ...
            'velocity_tolerance', velocity_tolerance, ...
            'acceleration_tolerance', acceleration_tolerance, ...
            'acceleration_assumption', 'instantaneous zero joint velocity', ...
            'current_velocity', [0; 0], 'damping', 0);
    end
end
