function [q, qdot, qddot, info] = coneik_inverse_state(gamma, gamma_dot, gamma_ddot, options)
%CONEIK_INVERSE_STATE Analytic position IK with checked second-order motion.
% All inputs/outputs use degrees, degrees/s and degrees/s^2, nominal joints.
% OPTIONS (optional scalar struct):
%   model           - coneik_model(alpha1_deg,alpha2_deg); default 20/20.
%   previous_q      - 2-vector; selects closest periodic analytic branch.
%                    At exact center, preserves this representative's q1.
%   current_velocity- 2-vector; an actual joint velocity, checked against
%                    gamma_dot. At singularities it preserves null motion.
%   on_infeasible   - 'error' (default) or 'report'. Report mode returns NaN
%                    for unsatisfied derivatives and info feasibility flags.
% Without options the original first analytic position branch is retained.
% Opt in with selection='constraint_aware' to evaluate both branches and
% bounded periodic representatives before choosing one. That mode also
% accepts joint_bounds, velocity_limits, acceleration_limits, center_q1,
% winding_radius and max_candidates; see coneik_select_inverse_state.
% In that selector, selection_policy='continuity_margin' adds soft
% preferred_velocity/preferred_acceleration, selection_weights [q v a margin]
% and selection_hysteresis. Defaults retain distance selection. Preferences
% never replace the hard measured current_velocity and do not prove path continuity.
% singular_motion='joint' (default) searches a rank-one velocity family when
% the minimum-norm state fails acceleration compatibility or motion limits.
% 'minimum_norm' retains the previous selected-velocity behavior. Optional
% velocity_limits/acceleration_limits apply to nominal joints in either mode.
% A provided current_velocity is never changed. The position representative
% remains selected by the rules above; this is not a trajectory/global proof.

    gamma = vector2(gamma, 'gamma');
    gamma_dot = vector2(gamma_dot, 'gamma_dot');
    gamma_ddot = vector2(gamma_ddot, 'gamma_ddot');
    if nargin < 4, options = struct(); end
    if ~(isstruct(options) && isscalar(options))
        error('coneik_inverse_state:InvalidOptions','options must be a scalar struct.');
    end
    if isfield(options,'selection')
        [q,qdot,qddot,info]=coneik_select_inverse_state(gamma,gamma_dot,gamma_ddot,options);
        return;
    end
    allowed = {'previous_q','current_velocity','on_infeasible','model',...
        'singular_motion','velocity_limits','acceleration_limits'};
    if any(~ismember(fieldnames(options), allowed))
        error('coneik_inverse_state:InvalidOptions','Unknown option field.');
    end
    options=coneik_motion_options(options);
    if isfield(options,'previous_q')
        if isempty(options.previous_q), options=rmfield(options,'previous_q');
        else, options.previous_q = vector2(options.previous_q,'previous_q'); end
    end
    if isfield(options,'current_velocity')
        if isempty(options.current_velocity), options=rmfield(options,'current_velocity');
        else, options.current_velocity = vector2(options.current_velocity,'current_velocity'); end
    end
    if ~isfield(options,'on_infeasible'), options.on_infeasible='error'; end
    policy=options.on_infeasible;
    if isstring(policy) && isscalar(policy), policy=char(policy); end
    if ~(ischar(policy) && isrow(policy) && any(strcmp(policy,{'error','report'})))
        error('coneik_inverse_state:InvalidOptions','on_infeasible must be error or report.');
    end

    model=[];
    if isfield(options,'model'), model=coneik_resolve_model(options.model); end
    [q11,q21,q12,q22,position_info] = inverse_motion_analysis(gamma(1),gamma(2),model);
    [q,qdot,qddot,info] = coneik_inverse_state_from_position(...
        gamma,gamma_dot,gamma_ddot,options,policy,model,[q11,q12;q21,q22],position_info);
end

function x=vector2(x,name)
    if ~(isnumeric(x) && isreal(x) && isvector(x) && numel(x)==2 && all(isfinite(x(:))))
        error('coneik_inverse_state:InvalidInput','%s must be a finite real numeric 2-vector.',name);
    end
    x=double(x(:));
end
