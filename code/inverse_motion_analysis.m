function [theta1_1, theta2_1, theta1_2, theta2_2, info] = inverse_motion_analysis(gamma1, gamma2, model)
%INVERSE_MOTION_ANALYSIS Stable inverse, optionally using CONEIK_MODEL geometry.
% An omitted/empty MODEL or coneik_model(20,20) selects the default geometry.
% Inputs: finite real scalar pitch/yaw in degrees (periodic angles accepted).
% Outputs: two NOMINAL angle pairs, in degrees, in the original branch order.
% The first four outputs preserve the existing calling convention. Optional
% INFO reports status, route, used_fallback, precision_dps and branch_count.
%
% Model: r = Rx(theta1)*Rz(alpha1)*Rx(theta2)*Rz(alpha2)*[1;0;0].
% theta2 is already nominal; do NOT subtract theta1 a second time.
% theta1 is in [-180,180). As in the original interface, theta1+theta2
% selects an actual-angle representative in [-180,180]; theta2 can exceed
% +/-180. Representatives at wrap cuts may differ by 360 degrees.
%
% Exact center: infinitely many solutions; retain (0,180),(0,-180).
% Exact outer boundary: the two outputs represent the same solution.
% Unreachable / invalid / unresolved inputs raise distinct errors.
% Ordinary cases and exact axial targets near the outer boundary use double
% precision. General difficult cases first try the optional native 50/80 then 100/200
% digit backend. If undecided/unavailable, the original Symbolic Math Toolbox
% ladder is retained. Both precision backends require two-level agreement.
% Guards are engineering dispatch rules, NOT certified error bounds.

    if nargin >= 3 && ~isempty(model)
        model = coneik_resolve_model(model);
        if ~model.is_default
            if nargout>4
                [theta1_1,theta2_1,theta1_2,theta2_2,info] = ...
                    coneik_inverse_position_resolved(gamma1,gamma2,model);
            else
                [theta1_1,theta2_1,theta1_2,theta2_2] = ...
                    coneik_inverse_position_resolved(gamma1,gamma2,model);
            end
            return;
        end
    end

    if nargout>4
        [theta1_1,theta2_1,theta1_2,theta2_2,info]=coneik_inverse_default_resolved(gamma1,gamma2,true);
    else
        [theta1_1,theta2_1,theta1_2,theta2_2]=coneik_inverse_default_resolved(gamma1,gamma2,true);
    end
end
