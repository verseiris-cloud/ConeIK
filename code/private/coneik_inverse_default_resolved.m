function [theta1_1,theta2_1,theta1_2,theta2_2,info]=coneik_inverse_default_resolved(gamma1,gamma2,allow_native_precision)
% Private default arithmetic; public validation/geometry dispatch stays outside.
    if ~(isnumeric(gamma1) && isreal(gamma1) && isscalar(gamma1) && ...
            isfinite(gamma1) && isnumeric(gamma2) && isreal(gamma2) && ...
            isscalar(gamma2) && isfinite(gamma2))
        error('inverse_motion_analysis:InvalidInput', ...
            'gamma1 and gamma2 must be finite real numeric scalars in degrees.');
    end
    g1 = double(gamma1); g2 = double(gamma2);
    if g1 < -180 || g1 >= 180, g1 = wrap_degrees(g1); end
    if g2 < -180 || g2 >= 180, g2 = wrap_degrees(g2); end
    if nargout > 4
        info = struct('status', 'regular', 'route', 'fast_half', ...
            'used_fallback', false, 'precision_dps', 16, 'branch_count', 2);
    end

    % Recognize exact degree identities, never a neighbourhood of a singularity.
    if (g1 == 0 && g2 == 0) || (g1 == -180 && g2 == -180)
        theta1_1 = 0; theta2_1 = 180;
        theta1_2 = 0; theta2_2 = -180;
        if nargout > 4
            info.status = 'center'; info.route = 'identity';
            info.branch_count = Inf;
        end
        return;
    end
    boundary_q1 = NaN;
    if g1 == 0 && abs(g2) == 40
        boundary_q1 = -180 * (g2 < 0);
    elseif g2 == 0 && abs(g1) == 40
        boundary_q1 = -90 * sign(g1);
    elseif g1 == -180 && abs(g2) == 140
        boundary_q1 = -180 * (g2 > 0);
    elseif g2 == -180 && abs(g1) == 140
        boundary_q1 = -90 * sign(g1);
    end
    if isfinite(boundary_q1)
        theta1_1 = boundary_q1; theta2_1 = 0;
        theta1_2 = boundary_q1; theta2_2 = 0;
        if nargout > 4
            info.status = 'outer_boundary'; info.route = 'identity';
            info.branch_count = 1;
        end
        return;
    end

    a = pi / 9; s = sin(a); c = cos(a);
    guard = 32 * eps; joint_budget_deg = 1e-8;
    if abs(g1) <= 90 && abs(g2) <= 90
        % Ordinary principal inputs need no degree argument reduction.
        p = g1 * (pi/180); y = g2 * (pi/180); cp = cos(p);
        x = cp*cos(y); ry = cp*sin(y); rz = -sin(p);
    else
        % Preserve quadrantal identities for nonprincipal periodic inputs.
        cp = cosd(g1);
        x = cp*cosd(g2); ry = cp*sind(g2); rz = -sind(g1);
    end
    rho = hypot(ry, rz);
    beta = atan2(rho, x);
    gap = 2*a - beta;
    product = sin(a + beta/2) * sin(a - beta/2);
    reason = '';
    if any(abs([g1,g2]) > 0 & abs([g1,g2]) < realmin*180/pi)
        reason = 'subnormal_input';
    elseif rho == 0 || any(abs([ry,rz,rho,beta]) > 0 & abs([ry,rz,rho,beta]) < realmin)
        reason = 'underflow';
    elseif abs(gap) <= guard
        reason = 'gap_guard';
    elseif gap < -guard
        error('inverse_motion_analysis:Unreachable', ...
            'Target direction is outside the fixed 20/20 degree reachable cone (beta > 40 degrees).');
    elseif ~isfinite(product) || product <= 0
        reason = 'radicand';
    elseif (180/pi)*(guard*cos(beta/2)/sqrt(product)) > joint_budget_deg
        reason = 'joint_budget';
    end

    if isempty(reason)
        u = sin(beta/2) / s;
        v = sqrt(product) / s;
        eta = (180/pi)*(2 * atan2(v, u));
        delta = (180/pi)*atan2(v, c*u);
        psi = (180/pi)*atan2(rz, ry);
        q = [psi+delta, -eta; psi-delta, eta];
    else
        native=false;
        [handled,q] = coneik_ikgeo_axis(g1,g2,20);
        if handled
            dps=16;reason='ikgeo_sp4_sp1_axis';
        else
            [q, dps, native] = precise_inverse(g1, g2, allow_native_precision);
        end
        if nargout > 4
            info.route = reason; if native,info.route='cpp_precision';end
            info.used_fallback = ~handled;
            info.precision_dps = dps;
        end
    end

    % Preserve the actual-angle convention without subtracting nearly equal
    % angles (which would erase a small nominal theta2 at the outer boundary).
    if ~all(isfinite(q(:)))
        error('inverse_motion_analysis:NumericalUnresolved', 'Nonfinite inverse output.');
    end
    theta1_1 = q(1,1); theta2_1 = q(1,2);
    theta1_2 = q(2,1); theta2_2 = q(2,2);
    % psi +/- delta lies within [-270,270], so one turn suffices. Scalar
    % branches avoid temporary logical-index arrays in the usual four-output API.
    if theta1_1 >= 180, theta1_1 = theta1_1 - 360;
    elseif theta1_1 < -180, theta1_1 = theta1_1 + 360; end
    if theta1_2 >= 180, theta1_2 = theta1_2 - 360;
    elseif theta1_2 < -180, theta1_2 = theta1_2 + 360; end
    if theta1_1 + theta2_1 > 180, theta2_1 = theta2_1 - 360;
    elseif theta1_1 + theta2_1 < -180, theta2_1 = theta2_1 + 360; end
    if theta1_2 + theta2_2 > 180, theta2_2 = theta2_2 - 360;
    elseif theta1_2 + theta2_2 < -180, theta2_2 = theta2_2 + 360; end
end

function [q, dps, native] = precise_inverse(g1, g2, allow_native_precision)
    native=false;
    if allow_native_precision
        [native,q,dps]=coneik_try_native_precision(g1,g2,[20,20]);
        if native,return;end
    end
    if exist('sym', 'file') == 0 || ~license('test', 'Symbolic_Toolbox')
        error('inverse_motion_analysis:PrecisionUnavailable', ...
            'This difficult input requires Symbolic Math Toolbox for reliable classification and angles.');
    end
    saved_digits = digits;
    restore_digits = onCleanup(@() digits(saved_digits)); %#ok<NASGU>
    % 'f' captures the exact binary64 rational, not a decimal approximation.
    g1_exact = sym(g1, 'f'); g2_exact = sym(g2, 'f');
    previous_state = ''; previous_q = [];
    for dps = [100, 200, 400, 800, 1024]
        digits(dps);
        rad = vpa(str2sym('pi')/180, dps);
        a = 20 * rad;
        p = vpa(g1_exact, dps) * rad; y = vpa(g2_exact, dps) * rad;
        x = cos(p)*cos(y); ry = cos(p)*sin(y); rz = -sin(p);
        beta = atan2(sqrt(ry*ry + rz*rz), x);
        gap = 2*a - beta;
        margin = vpa(sym(10)^(-dps+20), dps);
        if isAlways(abs(gap) <= margin) || isAlways(beta == 0)
            previous_state = ''; previous_q = [];
            continue;
        end
        if isAlways(gap < 0)
            state = 'unreachable'; q = [];
        else
            state = 'regular';
            s = sin(a); c = cos(a); h = beta/2;
            u = sin(h)/s;
            v = sqrt(sin(a+h)*sin(a-h))/s;
            eta = 2*atan2(v,u)/rad;
            delta = atan2(v,c*u)/rad;
            psi = atan2(rz,ry)/rad;
            q = double([psi+delta, -eta; psi-delta, eta]);
        end
        if strcmp(state, previous_state)
            if strcmp(state, 'unreachable')
                error('inverse_motion_analysis:Unreachable', ...
                    'Target direction is outside beta = 40 degrees (checked at increasing precision).');
            end
            difference = abs(wrap_degrees(q - previous_q));
            if all(isfinite(q(:))) && max(difference(:)) <= 1e-12
                return;
            end
        end
        previous_state = state; previous_q = q;
    end
    error('inverse_motion_analysis:NumericalUnresolved', ...
        'Classification or angles did not converge by 1024 decimal digits.');
end

function q = wrap_degrees(q)
    % Do not use mod(q+180,360)-180: the addition discards tiny angles.
    q = rem(q, 360);
    q(q >= 180) = q(q >= 180) - 360;
    q(q < -180) = q(q < -180) + 360;
    q(q == 0) = 0;
end
