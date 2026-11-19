function [q,info]=coneik_select_inverse_position_resolved(target,reference_q,c)
    model=c.model;bounds=c.bounds;weights=c.weights;weight_scale=c.weight_scale;w=c.w;
    objective=c.objective;policy=c.policy;tolerance=c.tolerance;
    detailed=nargout>1;
    if numel(target)==2
        [base,pos]=angle_position(target,model,tolerance,~detailed);
        g=rem(target,360);r=[cosd(g(1))*cosd(g(2)),cosd(g(1))*sind(g(2)),-sind(g(1))];
    else
        [base,pos]=coneik_direction_resolved(target,model,'report',tolerance,~detailed);
        r=target/max(abs(target));r=r/norm(r);
    end
    q=nan(2,1);
    candidate_angles=nan(2,2);candidate_feasible=false(1,2);candidate_cost=nan(1,2);
    selected_branch=NaN;normalized_objective=NaN;direction_error_deg=NaN;
    center=false;selection_complete=false;success=false;error_id='';error_message='';
    status='regular';
    if ~islogical(pos)
        status=pos.status{1};
        if ~pos.success
            error_id=pos.error_id{1};error_message=pos.error_message{1};
            finish();return;
        end
        center=strcmp(status,'center');
    end
    count=2;if center,count=1;end
    physical=inf(1,2);uncertain=false;
    for branch=1:count
        b=wrap(base(2*branch-1:2*branch)');candidate=nan(2,1);feasible=true;
        for axis=1:2
            if center&&axis==1
                [value,state]=free_angle(bounds(axis,:),reference_q(axis),objective,w(axis));
            else
                [value,state]=nearest_lift(b(axis),bounds(axis,:),reference_q(axis));
            end
            uncertain=uncertain||state<0;feasible=feasible&&state==1;candidate(axis)=value;
        end
        if ~feasible,continue;end
        phase_loss=abs(wrap(wrap(candidate)-b));if center,phase_loss(1)=0;end
        if any(phase_loss>min(1e-9,tolerance/8)),uncertain=true;continue;end
        delta=candidate-reference_q;
        physical(branch)=hypot(delta(1),delta(2));
        if strcmp(objective,'distance')
            score=hypot(sqrt(w(1))*delta(1),sqrt(w(2))*delta(2));
        else
            % Use original branch phases to avoid losing a small preference
            % difference merely through a large periodic representative.
            phase=b;if center,phase(1)=candidate(1);end
            d=wrap(wrap(phase)-wrap(reference_q));score=sum(w.*(4*sind(d/2).^2));
        end
        if ~isfinite(score)||~isfinite(physical(branch)),uncertain=true;continue;end
        candidate_angles(:,branch)=candidate;candidate_feasible(branch)=true;candidate_cost(branch)=score;
    end
    if uncertain
        status='numerical_unresolved';error_id='coneik_select_inverse_position:NumericalUnresolved';
        error_message='Periodic representatives or ranking exceed reliable numeric resolution.';
        finish();return;
    end
    feasible=find(candidate_feasible);
    if isempty(feasible)
        status='joint_bounds';selection_complete=true;
        error_id='coneik_select_inverse_position:NoBoundedSolution';error_message='No analytic branch has a representative inside joint_bounds.';
        finish();return;
    end
    best=feasible(1);
    if numel(feasible)==2
        second=feasible(2);
        if candidate_cost(second)<candidate_cost(best)||...
                (candidate_cost(second)==candidate_cost(best)&&physical(second)<physical(best))
            best=second;
        end
    end
    selected=candidate_angles(:,best);
    actual=forward(selected,model);residual=atan2d(norm(cross(actual,r)),dot(actual,r));
    direction_error_deg=residual;
    if ~isfinite(residual)||residual>tolerance
        status='residual_rejected';error_id='coneik_select_inverse_position:ResidualRejected';
        error_message='Selected representative exceeds direction_tolerance_deg.';
        finish();return;
    end
    q=selected;success=true;status='ok';selected_branch=best;
    normalized_objective=candidate_cost(best);selection_complete=true;
    finish();
    function finish()
        if detailed
            info=selection_report(c,reference_q,pos,base,candidate_angles,candidate_feasible,...
                candidate_cost,selected_branch,normalized_objective,direction_error_deg,...
                center,selection_complete,success,status,error_id,error_message);
        end
        if ~success && strcmp(policy,'error'),error(error_id,'%s',error_message);end
    end

end

function [Q,info]=angle_position(target,model,tolerance,compact)
% Reuse the installed binary64 kernel's existing acceptance guard. Only a
% confirmed ordinary success avoids batch-shaped diagnostic allocation;
% every other outcome runs the unchanged complete precision/failure route.
    backend=coneik_native_backend('position');
    if ~isempty(backend)
        [Q,code]=backend(target,model.alpha_deg);
        if code==0
            if compact,info=true;return;end
            row=struct('status','regular','route','cpp_analytic','used_fallback',false,...
                'precision_dps',16,'branch_count',2,'model_key',model.key);
            info=struct('model_key',model.key,'backend_requested','auto','backend_used','cpp',...
                'status',{{'regular'}},'success',true,'fallback',false,'used_precision',false,...
                'precision_known',true,'precision_dps',16,'error_id',{{''}},'error_message',{{''}},...
                'position_info',{{row}},'native_fast_count',1,'native_precision_count',0,'matlab_row_count',0);
            return;
        end
    end
    [Q,info]=coneik_inverse_batch_resolved(target,model,'auto','report','none',tolerance);
end

function [value,state]=nearest_lift(base,bounds,reference)
% State 1: accepted, 0: no bounded representative, -1: numeric ambiguity.
    value=NaN;state=0;
    low=-Inf;high=Inf;
    if isfinite(bounds(1)),low=ceil((bounds(1)-base)/360);end
    if isfinite(bounds(2)),high=floor((bounds(2)-base)/360);end
    % Expand only by adjacent integers to tolerate division rounding at an
    % endpoint; every returned double is checked against the actual bounds.
    near=round((reference-base)/360);near=min(max(near,low),high);
    % Do not substitute a distant, representable endpoint when the actual
    % closest turn cannot be ranked reliably in double precision.
    if ~isfinite(near)||abs(near)>2^40,state=-1;return;end
    turns=[near-1,near,near+1];
    turns=turns(abs(turns)<=2^40);
    if isempty(turns),state=-1;return;end
    values=base+360*turns;keep=values>=bounds(1)&values<=bounds(2);values=values(keep);
    if isempty(values),return;end
    % Turns are ascending: min's first-index tie break also chooses the
    % smaller physical angle without sorting or allocating all lifts.
    [~,best]=min(abs(values-reference));value=values(best);state=1;
end

function [value,state]=free_angle(bounds,reference,objective,weight)
    if strcmp(objective,'distance')||weight==0
        value=min(max(reference,bounds(1)),bounds(2));state=1;return;
    end
    [value,state]=nearest_lift(wrap(reference),bounds,reference);
    if state~=0,return;end
    ends=bounds(isfinite(bounds));
    if isempty(ends),state=-1;return;end
    differences=wrap(wrap(ends)-wrap(reference));cost=4*sind(differences/2).^2;
    [~,order]=sortrows([cost(:),abs(ends(:)-reference),ends(:)],[1,2,3]);
    value=ends(order(1));state=1;
end

function r=forward(q,m)
    q=rem(q,360);u=q(1);v=q(2);sa=m.sin_alpha(1);ca=m.cos_alpha(1);sb=m.sin_alpha(2);cb=m.cos_alpha(2);
    y=sa*cb+ca*sb*cosd(v);z=sb*sind(v);
    r=[ca*cb-sa*sb*cosd(v),cosd(u)*y-sind(u)*z,sind(u)*y+cosd(u)*z];
end
function x=wrap(x)
    x=rem(x,360);x(x>=180)=x(x>=180)-360;x(x< -180)=x(x< -180)+360;
end
function info=selection_report(c,reference_q,pos,base,angles,feasible,cost,best,score,residual,center,complete,success,status,error_id,error_message)
% Construct public diagnostics only when the caller requests them.
    info=struct('success',success,'status',status,'position_info',pos,'base_candidates',base,...
        'candidate_angles',angles,'candidate_feasible',feasible,'candidate_cost',cost,...
        'selected_branch',best,'objective',c.objective,'normalized_objective',score,'weights',c.weights,...
        'weight_scale',c.weight_scale,'reference_q',reference_q,'joint_bounds',c.bounds,...
        'direction_error_deg',residual,'direction_tolerance_deg',c.tolerance,'center_family_optimized',center,...
        'selection_complete',complete,'error_id',error_id,'error_message',error_message,...
        'scope','Separable static selection over analytic branches and periodic representatives; no motion or interval certificate.');
end
