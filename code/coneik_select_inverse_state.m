function [q,qdot,qddot,info] = coneik_select_inverse_state(gamma,gamma_dot,gamma_ddot,options)
%CONEIK_SELECT_INVERSE_STATE Evaluate analytic candidates before selecting one.
% Called by CONEIK_INVERSE_STATE with selection='constraint_aware'. Bounds are
% NOMINAL joint bounds (degree, degree/s, degree/s^2), not calibrated motor
% guarantees. joint_bounds is 2x2 [lower upper]; velocity_limits and
% acceleration_limits are nonnegative 2-vectors, allowing zero and Inf.
% Optional previous_q ranks passing candidates by physical angle distance.
% center_q1, winding_radius and max_candidates control position enumeration.
% The analytic base position solution is shared across representatives.
% Each representative retains the existing branch reconstruction and motion checks.
% Rank-one joint motion search is enabled by default; singular_motion can
% select the previous 'minimum_norm' behavior. No result from this finite
% position search establishes global or trajectory feasibility.
% selection_policy='continuity_margin' opts into normalized q/v/a preferences
% and nominal margin scoring. preferred_velocity/preferred_acceleration are
% soft finite 2-vectors, distinct from hard current_velocity. selection_weights
% is [position velocity acceleration margin] (default ones); selection_hysteresis
% is a nonnegative score deadband retaining the nearest feasible q. Omitted
% policy or 'distance' preserves the original outputs. See task usage notes.

    if nargin<4, options=struct(); end
    allowed={'selection','previous_q','current_velocity','on_infeasible','joint_bounds',...
        'velocity_limits','acceleration_limits','center_q1','max_candidates','winding_radius','model','singular_motion',...
        'selection_policy','preferred_velocity','preferred_acceleration','selection_weights','selection_hysteresis'};
    if ~(isstruct(options)&&isscalar(options)) || any(~ismember(fieldnames(options),allowed))
        error('coneik_select_inverse_state:InvalidOptions','Unknown option or non-scalar options struct.');
    end
    if isfield(options,'selection') && ~is_text(options.selection,'constraint_aware')
        error('coneik_select_inverse_state:InvalidOptions','selection must be constraint_aware.');
    end
    options=coneik_motion_options(options);
    continuity=false;
    if isfield(options,'selection_policy')
        if is_text(options.selection_policy,'continuity_margin'),continuity=true;
        elseif ~is_text(options.selection_policy,'distance')
            error('coneik_select_inverse_state:InvalidOptions','selection_policy must be distance or continuity_margin.');
        end
    end
    tuning={'preferred_velocity','preferred_acceleration','selection_weights','selection_hysteresis'};
    if ~continuity&&any(isfield(options,tuning))
        error('coneik_select_inverse_state:InvalidOptions','Soft selection options require continuity_margin.');
    end
    policy='error';
    if isfield(options,'on_infeasible'), policy=options.on_infeasible; end
    if ~(is_text(policy,'error')||is_text(policy,'report'))
        error('coneik_select_inverse_state:InvalidOptions','on_infeasible must be error or report.');
    end
    velocity_limits=limit_option(options,'velocity_limits');
    acceleration_limits=limit_option(options,'acceleration_limits');
    previous=[];
    if isfield(options,'previous_q') && ~isempty(options.previous_q)
        previous=vector2(options.previous_q,'previous_q');
    end
    current_velocity=[];
    if isfield(options,'current_velocity') && ~isempty(options.current_velocity)
        current_velocity=vector2(options.current_velocity,'current_velocity');
    end
    gamma=vector2(gamma,'gamma'); gamma_dot=vector2(gamma_dot,'gamma_dot'); gamma_ddot=vector2(gamma_ddot,'gamma_ddot');
    position_options=struct();
    fields={'previous_q','joint_bounds','center_q1','max_candidates','winding_radius','model'};
    for k=1:numel(fields)
        if isfield(options,fields{k}), position_options.(fields{k})=options.(fields{k}); end
    end
    % Preserve the public candidate adapter's validation order, resolving
    % this call's explicit model once before enumeration-specific checks.
    model=[];
    if isfield(options,'model'),model=coneik_resolve_model(options.model);end
    [Q,search]=coneik_position_candidates_resolved(gamma,position_options,[],model);
    tolerances=struct('position',1e-9,'velocity',1e-8,'acceleration',1e-7);
    bounds=search.joint_bounds;
    if continuity,context=selection_context(options,previous,bounds,velocity_limits,acceleration_limits);end
    template=struct('q',nan(2,1),'qdot',nan(2,1),'qddot',nan(2,1),'branch',NaN,...
        'accepted',false,'reason','','rank',NaN,'motion_feasible',false,...
        'position_residual',NaN,'velocity_residual',NaN,'acceleration_residual',NaN,...
        'velocity_source','',...
        'position_violation',nan(2,1),'velocity_violation',nan(2,1),...
        'acceleration_violation',nan(2,1),'selection_cost',Inf);
    reports=repmat(template,1,size(Q,2)); states=cell(1,size(Q,2));
    if continuity
        for k=1:numel(reports),reports(k).selection_score=Inf;reports(k).score_detail=[];end
    end
    % Avoid cache overhead when the enumerated set has little repetition.
    % This is only an enablement heuristic: actual hits still use the legacy
    % reconstructed q inside the evaluator, with exact signed-zero bits.
    reuse=[];
    if size(Q,2)>2
        keys=reshape(typecast(full(rem(Q(:),360)),'uint64'),2,[]);
        if size(Q,2)<=32
            % Same unique-key count, without sorting small candidate sets.
            % Column j is a duplicate iff an earlier column has both bits.
            duplicate=any(triu(keys(1,:)'==keys(1,:) & keys(2,:)'==keys(2,:),1),1);
            unique_count=nnz(~duplicate);
        else
            unique_count=size(unique(keys','rows'),1);
        end
        if unique_count<=size(Q,2)/2
            reuse=struct('keys',zeros(2,0,'uint64'),'geometry',{{}});
        end
    end
    % Only previous_q and continuity fixed_q vary between candidates.
    % Preserve option field order while preparing shared options once.
    legacy=struct('previous_q',[],'on_infeasible','report');
    legacy.velocity_limits=velocity_limits;legacy.acceleration_limits=acceleration_limits;
    legacy.joint_search_limits_only=true;
    if isfield(options,'singular_motion'),legacy.singular_motion=options.singular_motion;end
    if isfield(options,'model'),legacy.model=options.model;end
    if ~isempty(current_velocity),legacy.current_velocity=current_velocity;end
    if continuity,legacy.state_selection_context=context;end
    for k=1:size(Q,2)
        reports(k).q=Q(:,k); reports(k).branch=search.branch(k);
        if isempty(previous), reports(k).selection_cost=k;
        else, reports(k).selection_cost=norm(Q(:,k)-previous); end
        legacy.previous_q=Q(:,k);
        if continuity
            legacy.state_selection_context.fixed_q=Q(:,k);
        end
        try
            [selected,v,a,state,reuse]=coneik_inverse_state_resolved(...
                gamma,gamma_dot,gamma_ddot,legacy,'report',model,...
                search.base_candidates,search.position_info,reuse);
        catch exception
            if ~strcmp(exception.identifier,'coneik_inverse_state:NumericalUnresolved'), rethrow(exception); end
            reports(k).reason='numerical_unresolved';
            continue;
        end
        if norm(selected-Q(:,k),Inf)>tolerances.position
            error('coneik_select_inverse_state:CandidateMismatch','Legacy solver did not retain the intended physical representative.');
        end
        states{k}=state; reports(k).qdot=v; reports(k).qddot=a;
        reports(k).rank=state.rank; reports(k).motion_feasible=state.motion_feasible;
        reports(k).position_residual=state.position_residual;
        reports(k).velocity_residual=state.velocity_residual;
        reports(k).acceleration_residual=state.acceleration_residual;
        reports(k).velocity_source=state.velocity_source;
        reports(k).position_violation=max([bounds(:,1)-selected,selected-bounds(:,2),zeros(2,1)],[],2);
        reports(k).velocity_violation=max(abs(v)-velocity_limits,0);
        reports(k).acceleration_violation=max(abs(a)-acceleration_limits,0);
        if ~state.motion_feasible, reports(k).reason=state.failure_reason;
        elseif any(reports(k).position_violation>tolerances.position), reports(k).reason='joint_bounds';
        elseif any(reports(k).velocity_violation>tolerances.velocity), reports(k).reason='velocity_limits';
        elseif any(reports(k).acceleration_violation>tolerances.acceleration), reports(k).reason='acceleration_limits';
        else, reports(k).accepted=true; reports(k).reason='accepted'; end
        if continuity&&reports(k).accepted
            [reports(k).selection_score,reports(k).score_detail]=coneik_state_selection_score(context,selected,v,a);
        end
    end
    passing=find([reports.accepted]); selected_index=0;
    best_index=0;incumbent=0;retained=false;
    has_unresolved=any(strcmp({reports.reason},'numerical_unresolved'));
    if isempty(passing)
        q=nan(2,1); qdot=q; qddot=q;
        info=struct('status','no_feasible_candidate','position_info',search.position_info,...
            'branch',NaN,'selected_q',q,'velocity_feasible',false,'acceleration_feasible',false,...
            'motion_feasible',false,'failure_reason','NoFeasibleCandidate');
        if has_unresolved
            info.status='numerical_unresolved';
            info.failure_reason='NumericalUnresolvedCandidates';
        end
    else
        [~,best]=min([reports(passing).selection_cost]); selected_index=passing(best);
        if continuity
            incumbent=selected_index;
            table=[[reports(passing).selection_score]',[reports(passing).selection_cost]',passing'];
            [~,order]=sortrows(table,[1 2 3]);best_index=passing(order(1));
            selected_index=best_index;
            if ~isempty(previous)&&reports(incumbent).selection_score-reports(best_index).selection_score<=context.hysteresis
                selected_index=incumbent;retained=incumbent~=best_index;
            end
        end
        q=reports(selected_index).q; qdot=reports(selected_index).qdot; qddot=reports(selected_index).qddot;
        info=states{selected_index}; info.branch=reports(selected_index).branch;
    end
    info.selection='constraint_aware'; info.candidates=reports; info.candidate_search=search;
    if continuity
        info.selection_policy='continuity_margin';
        info.selection_context=context;
        info.selection_ranking=struct('best_candidate',best_index,'incumbent_candidate',incumbent,...
            'hysteresis_retained',retained,'selected_score',NaN);
        if selected_index>0,info.selection_ranking.selected_score=reports(selected_index).selection_score;end
    end
    if isfield(search,'model_key'), info.model_key=search.model_key; end
    info.candidate_count=numel(reports); info.accepted_candidate_count=numel(passing);
    info.selected_candidate=selected_index; info.constraints_feasible=~isempty(passing);
    info.joint_bounds=bounds; info.velocity_limits=velocity_limits; info.acceleration_limits=acceleration_limits;
    info.constraint_tolerances=tolerances;
    info.search_complete=search.position_enumeration_complete && ~has_unresolved && all([reports.rank]==2);
    info.singular_nullspaces_exhausted=false; info.global_feasibility_checked=false;
    if any(cellfun(@(s)isstruct(s)&&isfield(s,'joint_motion'),states))
        info.joint_motion_searches=cell(size(states));
        for k=1:numel(states)
            if isstruct(states{k})&&isfield(states{k},'joint_motion')
                info.joint_motion_searches{k}=states{k}.joint_motion;
            end
        end
    end
    info.feasibility_scope='tested nominal candidates and selected singular velocity/acceleration; no trajectory or global proof';
    if isfield(info,'joint_motion_searches')
        info.feasibility_scope='tested nominal candidates with rank-one joint motion recovery; no trajectory or global proof';
    end
    if isempty(passing) && is_text(policy,'error')
        if has_unresolved
            error('coneik_inverse_state:NumericalUnresolvedCandidates',...
                'No passing candidate was found and at least one candidate could not be resolved numerically. This does not establish infeasibility.');
        end
        error('coneik_inverse_state:NoFeasibleCandidate',...
            'No tested analytic candidate satisfies the motion and nominal bounds. This finite search is not a global infeasibility proof.');
    end
end

function context=selection_context(options,previous,bounds,vl,al)
    context=struct('previous_q',previous,'preferred_velocity',[],...
        'preferred_acceleration',[],'bounds',bounds,'velocity_limits',vl,...
        'acceleration_limits',al,'weights',ones(1,4),'hysteresis',0);
    for name={'preferred_velocity','preferred_acceleration'}
        key=name{1};if ~isfield(options,key),continue;end
        value=options.(key);
        if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&all(isfinite(value(:))))
            error('coneik_select_inverse_state:InvalidOptions','%s must be a finite real 2-vector.',key);
        end
        context.(key)=double(value(:));
    end
    if isfield(options,'selection_weights')
        w=options.selection_weights;
        if ~(isnumeric(w)&&isreal(w)&&isvector(w)&&numel(w)==4&&all(isfinite(w(:)))&&all(w(:)>=0)&&any(w(:)>0))
            error('coneik_select_inverse_state:InvalidOptions','selection_weights must be a finite nonnegative nonzero 4-vector.');
        end
        context.weights=double(w(:)');
    end
    if isfield(options,'selection_hysteresis')
        h=options.selection_hysteresis;
        if ~(isnumeric(h)&&isreal(h)&&isscalar(h)&&isfinite(h)&&h>=0)
            error('coneik_select_inverse_state:InvalidOptions','selection_hysteresis must be finite and nonnegative.');
        end
        context.hysteresis=double(h);
    end
    span=bounds(:,2)-bounds(:,1);context.position_scale=360*ones(2,1);
    finite_span=isfinite(span)&span>0;context.position_scale(finite_span)=span(finite_span);
    context.velocity_scale=motion_scale(vl,context.preferred_velocity);
    context.acceleration_scale=motion_scale(al,context.preferred_acceleration);
end

function scale=motion_scale(limit,preference)
    scale=ones(2,1);if ~isempty(preference),scale=max(scale,abs(preference));end
    finite_limit=isfinite(limit)&limit>0;scale(finite_limit)=limit(finite_limit);
end

function yes=is_text(value,expected)
    yes=(ischar(value)&&isrow(value)&&strcmp(value,expected))||...
        (isstring(value)&&isscalar(value)&&~ismissing(value)&&strcmp(value,expected));
end

function value=limit_option(options,name)
    value=[Inf;Inf];
    if isfield(options,name), value=options.(name); end
    if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&~any(isnan(value(:)))&&all(value(:)>=0))
        error('coneik_select_inverse_state:InvalidOptions','%s must be a nonnegative real 2-vector; Inf is permitted.',name);
    end
    value=double(value(:));
end

function value=vector2(value,name)
    if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&all(isfinite(value(:))))
        error('coneik_select_inverse_state:InvalidInput','%s must be a finite real 2-vector.',name);
    end
    value=double(value(:));
end
