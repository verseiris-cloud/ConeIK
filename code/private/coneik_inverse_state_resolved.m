function [q,qdot,qddot,info,reuse] = coneik_inverse_state_resolved(gamma,gamma_dot,gamma_ddot,options,policy,model,candidates,position_info,reuse)
% Private post-IK evaluator. The caller validates inputs/model first.
% REUSE belongs to exactly one selector invocation with a fixed model.
% Only geometry and its SVD are reused; every candidate motion is evaluated.
    if nargin<9,reuse=[];end
    if cosd(rem(gamma(1),360))==0
        error('coneik_inverse_state:ChartSingularity','Pitch/yaw motion is undefined at an exact direction pole.');
    end
    branch=1;
    q = candidates(:,1);
    if isfield(options,'previous_q')
        previous = options.previous_q;
        if strcmp(position_info.status,'center')
            q = [previous(1);180+360*round((previous(2)-180)/360)];
            branch=0; % representative of a continuous family, not a branch
        else
            lifted = candidates + 360*round((previous-candidates)/360);
            [~,branch] = min(sum((lifted-previous).^2,1));
            q = lifted(:,branch);
        end
    end
    if ~all(isfinite(q))
        error('coneik_inverse_state:NumericalUnresolved','Nonfinite periodic representative.');
    end

    % The position IK accepts nonprincipal pitch charts. For cos(pitch)<0,
    % canonical pitch equals +/-180-pitch: both its derivatives change sign.
    chart_sign = 1;
    if cosd(rem(gamma(1),360)) < 0, chart_sign = -1; end
    gd = gamma_dot; ga = gamma_ddot;
    gd(1)=chart_sign*gd(1); ga(1)=chart_sign*ga(1);
    hit=[];
    if ~isempty(reuse)
        % Key the actual legacy-reconstructed q, not the enumerated Q.
        % Integer bit comparisons keep signed zero and rounded lifts apart.
        key=typecast(full(rem(double(q(:)),360)),'uint64');
        hit=find(all(reuse.keys==key,1),1);
    end
    if isempty(hit)
        [r,canonical_gamma,J,~,detail]=coneik_kinematics_resolved(q,[0;0],model);
        factor=coneik_jacobian_state(J);
        if ~isempty(reuse)
            reuse.keys(:,end+1)=key;
            reuse.geometry{end+1}=struct('r',r,'canonical_gamma',canonical_gamma,...
                'J',J,'detail',detail,'factor',factor);
        end
    else
        geometry=reuse.geometry{hit};
        r=geometry.r;canonical_gamma=geometry.canonical_gamma;J=geometry.J;
        detail=geometry.detail;factor=geometry.factor;
    end
    target=[cosd(gamma(1))*cosd(gamma(2));cosd(gamma(1))*sind(gamma(2));-sind(gamma(1))];
    position_residual=norm(r-target);
    if ~isfinite(position_residual) || position_residual>1e-10
        error('coneik_inverse_state:NumericalUnresolved','Position residual exceeds 1e-10 in unit-direction chord length.');
    end
    status='regular'; solver='backslash';
    if factor.rank<2, status='singular_nonunique'; solver='retained_svd'; end
    info=struct('status',status,'position_info',position_info,'branch',branch,...
        'selected_q',q,'canonical_gamma',canonical_gamma,'pitch_chart_sign',chart_sign,...
        'task_chart','principal gamma','position_residual',position_residual,...
        'position_tolerance',1e-10,'rank',factor.rank,'singular_values',factor.s,...
        'condition',factor.condition,'near_singular',factor.near_singular,...
        'solver',solver,'damping',0,'velocity_source','',...
        'target_velocity',gd,'target_acceleration',ga,...
        'achieved_velocity',nan(2,1),'achieved_acceleration',nan(2,1),...
        'curvature',nan(2,1),'velocity_residual',NaN,'acceleration_residual',NaN,...
        'velocity_tolerance',1e-8+1e-10*norm(gd,Inf),...
        'acceleration_tolerance',1e-7+1e-10*norm(ga,Inf),...
        'velocity_feasible',false,'acceleration_feasible',false,'motion_feasible',false,...
        'feasibility_scope','selected configuration and velocity',...
        'global_feasibility_checked',false,'failure_reason','');
    if isfield(options,'model'), info.model_key=model.key; end
    if isfield(options,'current_velocity')
        qdot=options.current_velocity;
        velocity_source='provided';
    else
        qdot=solve(J,gd,factor);
        velocity_source='minimum_norm';
        if factor.rank==2, velocity_source='unique'; end
    end
    info.velocity_source=velocity_source;
    info.achieved_velocity=J*qdot;
    info.velocity_residual=norm(info.achieved_velocity-gd,Inf);
    if ~all(isfinite(qdot)) || ~isfinite(info.velocity_residual)
        error('coneik_inverse_state:NumericalUnresolved','Nonfinite velocity calculation.');
    end
    if info.velocity_residual>info.velocity_tolerance
        if factor.rank==2 && ~isfield(options,'current_velocity')
            error('coneik_inverse_state:NumericalUnresolved','Full-rank velocity solve did not meet the task residual tolerance.');
        end
        info.status='velocity_infeasible'; info.failure_reason='VelocityInfeasibleAtSelectedConfiguration';
        if strcmp(policy,'error')
            error('coneik_inverse_state:VelocityInfeasibleAtSelectedConfiguration',...
                'Requested velocity is not met at the selected configuration (residual %.6g deg/s).',info.velocity_residual);
        end
        qdot=nan(2,1); qddot=nan(2,1); return;
    end
    info.velocity_feasible=true;
    curvature=detail.curvature_coefficients*[qdot(1)^2;qdot(1)*qdot(2);qdot(2)^2];
    rhs=ga-curvature;
    qddot=solve(J,rhs,factor);
    velocity_limits=[Inf;Inf];acceleration_limits=[Inf;Inf];
    if isfield(options,'velocity_limits'),velocity_limits=options.velocity_limits;end
    if isfield(options,'acceleration_limits'),acceleration_limits=options.acceleration_limits;end
    joint_mode=~isfield(options,'singular_motion')||strcmp(options.singular_motion,'joint');
    finite_motion=all(isfinite([curvature;qddot]));
    acceleration_ok=finite_motion && norm(J*qddot+curvature-ga,Inf)<=info.acceleration_tolerance;
    limits_ok=all(abs(qdot)<=velocity_limits+1e-8)&&all(abs(qddot)<=acceleration_limits+1e-7);
    refining=isfield(options,'state_selection_context');
    baseline_feasible=acceleration_ok&&limits_ok;
    if factor.rank==1 && joint_mode && finite_motion && (~acceleration_ok||~limits_ok||refining)
        try
            [joint_v,joint_a,joint_report]=coneik_singular_motion(J,detail.curvature_coefficients,...
                gd,ga,factor,qdot,options,info.velocity_tolerance,info.acceleration_tolerance);
        catch ex
            if ~(refining&&baseline_feasible&&strcmp(ex.identifier,'coneik_inverse_state:NumericalUnresolved')),rethrow(ex);end
            joint_report=struct('attempted',true,'status','optional_refinement_unresolved',...
                'failure_identifier',ex.identifier,'certified',false,'global_feasibility_checked',false);
            joint_v=qdot;joint_a=qddot;
        end
        fallback=false;
        if refining&&baseline_feasible
            baseline_score=coneik_state_selection_score(options.state_selection_context,[],qdot,qddot);
            fallback=~strcmp(joint_report.status,'feasible');
            if ~fallback
                refined_score=coneik_state_selection_score(options.state_selection_context,[],joint_v,joint_a);
                fallback=baseline_score<=refined_score;
            end
            joint_report.original_feasible_fallback_retained=fallback;
            joint_report.original_feasible_score=baseline_score;
        end
        info.joint_motion=joint_report;
        info.feasibility_scope='rank-one motion family at selected configuration; numerical, not global or trajectory proof';
        if fallback
            % Optional soft refinement cannot invalidate an existing witness.
        elseif strcmp(joint_report.status,'feasible')
            qdot=joint_v;qddot=joint_a;
            if ~isfield(options,'current_velocity'),info.velocity_source='joint_nullspace';end
            info.solver='retained_svd_joint';
            info.achieved_velocity=J*qdot;
            info.velocity_residual=norm(info.achieved_velocity-gd,Inf);
            curvature=detail.curvature_coefficients*[qdot(1)^2;qdot(1)*qdot(2);qdot(2)^2];
            limits_ok=true;
        else
            info.status='joint_motion_unavailable';info.failure_reason='NoJointMotionAtSelectedConfiguration';
            info.curvature=curvature;info.acceleration_residual=norm(J*qddot+curvature-ga,Inf);
            if strcmp(policy,'error')
                error('coneik_inverse_state:NoJointMotionAtSelectedConfiguration',...
                    'No joint motion witness was found at the selected configuration. This numerical search is not a global infeasibility proof.');
            end
            qddot=nan(2,1);return;
        end
    end
    info.curvature=curvature;
    info.achieved_acceleration=J*qddot+curvature;
    info.acceleration_residual=norm(info.achieved_acceleration-ga,Inf);
    if ~all(isfinite([curvature;qddot])) || ~isfinite(info.acceleration_residual)
        error('coneik_inverse_state:NumericalUnresolved','Nonfinite acceleration calculation.');
    end
    if info.acceleration_residual>info.acceleration_tolerance
        if factor.rank==2
            error('coneik_inverse_state:NumericalUnresolved','Full-rank acceleration solve did not meet the task residual tolerance.');
        end
        info.status='acceleration_infeasible'; info.failure_reason='AccelerationInfeasibleAtSelectedVelocity';
        if strcmp(policy,'error')
            error('coneik_inverse_state:AccelerationInfeasibleAtSelectedVelocity',...
                'Requested acceleration is not met for the selected joint velocity (residual %.6g deg/s^2). No null-velocity search was performed.',info.acceleration_residual);
        end
        qddot=nan(2,1); return;
    end
    info.acceleration_feasible=true;
    if ~limits_ok && ~isfield(options,'joint_search_limits_only')
        info.status='motion_limits_infeasible';info.failure_reason='MotionLimitsAtSelectedState';
        if strcmp(policy,'error')
            error('coneik_inverse_state:MotionLimitsAtSelectedState','Selected joint state violates the supplied motion limits.');
        end
        return;
    end
    info.motion_feasible=true;
end

function x=solve(J,b,f)
    if f.rank==2, x=J\b;
    elseif any(f.active), x=f.V(:,f.active)*((f.U(:,f.active)'*b)./f.s(f.active));
    else, x=zeros(2,1); end
end

