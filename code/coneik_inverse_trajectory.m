function [result, stats, trajectory] = coneik_inverse_trajectory(t, gamma, options)
%CONEIK_INVERSE_TRAJECTORY Finite candidate beam search followed by cubic checks.
% [result,stats,trajectory]=coneik_inverse_trajectory(t,gamma,options).
% T is a finite strictly increasing vector of seconds; GAMMA is N-by-2
% pitch/yaw in degrees. Nonfinite/unreachable targets split segments.
% Result columns: [t,q1,q2,v1,v2,a1,a2], all nominal degree coordinates.
% Each retained full path is interpolated by one not-a-knot spline per axis
% (two points: linear; three: quadratic) and checked over every interval.
% The same coefficients provide all six output columns and later evaluation.
% OPTIONS: mode='continuous', joint_bounds (2x2 [lower upper]),
% velocity_limits, acceleration_limits (nonnegative 2-vectors), previous_q
% (preference only), beam_width (8), max_candidates (256), winding_radius (1),
% center_q1 (finite representative seeds). Defaults impose no joint bounds.
% model=coneik_model(alpha1,alpha2) sets known cone angles in degrees. Omit it
% for the original 20/20 model and the original output structure.
% position_backend='matlab'/'cpp'/'auto' opts into one batch position call.
% It only changes cached position geometry, not beam search or verification.
% Omit position_backend to preserve the original per-row MATLAB route.
% Candidates are selected in displacement order after hard verification.
% full_candidate_diagnostics=false stops at the first fully verified path.
% Set true to inspect all retained paths; this does not change selection.
% not_checked entries have empty reports; passing counts cover checked paths.
% Initial velocity is free. Singleton segments fail. No gap is bridged.
% This finite search is incomplete: failure is not global infeasibility.
% Optional intersample_direction_tolerance_deg (0<degrees<180) checks the
% shortest-arc SLERP of target directions. The optional integer
% intersample_direction_max_refinements defaults to 1024/segment. Bisections
% verify the stored curve; they do not add waypoints or alter its timing.
% Optional initial_q fixes the physical lifted start at source t(1).
% initial_velocity / final_velocity fix derivatives at source endpoints.
% joint_excursion_bounds is a finite 2x2 signed offset window from initial_q
% (lower<=0<=upper), not accumulated travel. It applies across stored segments;
% endpoint requests are never relocated across gaps. Endpoint acceleration is free.
    if nargin<3 || isempty(options),options=struct();end
    opts=validate_options(options);
    resolved_model=[];if isfield(opts,'model'),resolved_model=opts.model;end
    if ~(isnumeric(t)&&isreal(t)&&isvector(t)&&~isempty(t)&&all(isfinite(t(:))))
        error('coneik_inverse_trajectory:InvalidTime','t must be a nonempty finite real vector.');
    end
    t=double(t(:));n=numel(t);
    if any(~isfinite(diff(t))) || any(diff(t)<=0)
        error('coneik_inverse_trajectory:InvalidTime','Time must be strictly increasing.');
    end
    if ~(isnumeric(gamma)&&isreal(gamma)&&isequal(size(gamma),[n 2]))
        error('coneik_inverse_trajectory:InvalidGamma','gamma must be N-by-2 real numeric data.');
    end
    gamma=double(gamma);
    result=[t,nan(n,6)];cache=cell(n,1);valid=false(n,1);
    row_status=repmat({'target_gap'},n,1);row_error=repmat({''},n,1);
    cache_options=struct('winding_radius',0,'max_candidates',2);
    if isfield(opts,'model'),cache_options.model=opts.model;end
    % Cache unconstrained position geometry once per sample. Bounds are not
    % used here: a constrained search failure must not subdivide a segment.
    if isfield(opts,'position_backend')
        batch_model=[];if isfield(opts,'model'),batch_model=opts.model;end
        [batch_q,batch_info]=coneik_inverse_batch(gamma,batch_model,...
            struct('backend',opts.position_backend,'on_error','report'));
        for i=1:n
            if any(~isfinite(gamma(i,:)))
                row_error{i}='nonfinite_target';continue;
            end
            if batch_info.success(i)
                cache{i}=struct('input_gamma',gamma(i,:)',...
                    'base_candidates',reshape(batch_q(i,:),2,2),...
                    'position_info',batch_info.position_info{i});
                if isfield(opts,'model'),cache{i}.model_key=opts.model.key;end
                valid(i)=true;row_status{i}='pending_search';
            else
                row_error{i}=batch_info.error_id{i};
                if strcmp(row_error{i},'inverse_motion_analysis:Unreachable')
                    row_status{i}='unreachable_target';
                else
                    row_status{i}='target_unresolved';
                end
            end
        end
    else
        for i=1:n
            if any(~isfinite(gamma(i,:)))
                row_error{i}='nonfinite_target';continue;
            end
            try
                [~,meta]=coneik_position_candidates_resolved(gamma(i,:)',cache_options,[],resolved_model);
                cache{i}=meta.cache;valid(i)=true;row_status{i}='pending_search';
            catch exception
                row_error{i}=exception.identifier;
                if strcmp(exception.identifier,'inverse_motion_analysis:Unreachable')
                    row_status{i}='unreachable_target';
                else
                    row_status{i}='target_unresolved';
                end
            end
        end
    end
    starts=find(diff([false;valid])==1);ends=find(diff([valid;false])==-1);
    trajectory=struct('schema_version',1,'mode','continuous',...
        'coordinate_system','nominal degrees','time_unit','seconds',...
        'interpolation','not-a-knot cubic; two points linear; three points quadratic',...
        'initial_velocity','free','previous_q_role','preference only, not a boundary condition',...
        'intersample_task_tracking_verified',false,'global_feasibility_checked',false,...
        'segments',struct([]));
    if isfield(opts,'model')
        trajectory.schema_version=2;
        trajectory.model=opts.model;trajectory.model_key=opts.model.key;
    end
    segments=struct([]);
    for k=1:numel(starts)
        idx=starts(k):ends(k);
        [segment,entry]=solve_segment(t(idx),gamma(idx,:),cache(idx),idx,opts,n);
        if k==1,segments=entry;else,segments(k)=entry;end
        if strcmp(entry.status,'pass')
            if isempty(trajectory.segments),trajectory.segments=segment;
            else,trajectory.segments(end+1)=segment;end %#ok<AGROW>
            part=struct('segments',segment);
            result(idx,:)=coneik_trajectory_eval(part,t(idx));
            row_status(idx)={'pass'};
        else
            row_status(idx)={entry.status};
            row_error(idx)={entry.failure_reason};
        end
    end
    failed=any(~isfinite(result(:,2:7)),2);
    stats=struct('n_total',n,'n_ok',nnz(~failed),'n_fail',nnz(failed),...
        'fail_rows',find(failed)','mode','continuous','row_status',{row_status},...
        'row_error',{row_error},'segments',segments,...
        'segment_count',numel(segments),'accepted_segment_count',numel(trajectory.segments),...
        'target_fail_rows',find(~valid)','search_fail_rows',find(valid&failed)',...
        'motion_feasible',~any(failed),'global_feasibility_checked',false,...
        'search_complete',false,'beam_width',opts.beam_width,'max_candidates',opts.max_candidates,...
        'winding_radius',opts.winding_radius,'initial_velocity','free',...
        'previous_q_role','preference only, not a boundary condition',...
        'intersample_task_tracking_verified',false,...
        'verification_method','float64 polynomial endpoints and stationary roots',...
        'certified',false,'status','pass');
    if any(failed)
        if stats.n_ok>0,stats.status='partial';else,stats.status='no_accepted_trajectory';end
    end
    stats.continuous_joint_bounds_verified=~any(failed);
    trajectory.source_t=t;trajectory.options=opts;
    if isfield(opts,'model'),stats.model_key=opts.model.key;end
    if isfield(opts,'position_backend')
        stats.position_backend=opts.position_backend;
        stats.position_batch_summary=struct('backend_used',batch_info.backend_used,...
            'rows',n,'success_rows',nnz(batch_info.success),'failed_rows',nnz(~batch_info.success),...
            'native_fast_count',batch_info.native_fast_count,'matlab_row_count',batch_info.matlab_row_count,...
            'fallback_rows',nnz(batch_info.fallback),'precision_rows',nnz(batch_info.used_precision),...
            'scope','position geometry only; beam search and verification unchanged');
    end
    stats.verification=coneik_verify_trajectory(trajectory,opts,result);
    boundary_enabled=any(isfield(opts,{'initial_q','initial_velocity','final_velocity','joint_excursion_bounds'}));
    coverage_only=false;
    if boundary_enabled
        stats.boundary_constraints=stats.verification.boundary_constraints;
        stats.boundary_constraints_verified=stats.boundary_constraints.verified;
        trajectory.boundary_constraints_verified=stats.boundary_constraints_verified;
        stats.boundary_scope='hard q/v at original source endpoints; signed excursion is an absolute physical window, not accumulated travel; gaps not bridged';
        trajectory.boundary_scope=stats.boundary_scope;
        coverage_only=stats.boundary_constraints.coverage_only_failure;
        if isfield(opts,'initial_velocity')
            stats.initial_velocity='hard_requested';trajectory.initial_velocity='hard_requested';
        end
        if isfield(opts,'initial_velocity')||isfield(opts,'final_velocity')
            trajectory.interpolation='segment-specific C2 cubic; requested source-endpoint velocities; otherwise legacy not-a-knot';
        end
    end
    if isfield(opts,'intersample_direction_tolerance_deg')
        stats.intersample_task_tracking_verified=stats.verification.intersample_task_tracking_verified;
        trajectory.intersample_task_tracking_verified=stats.intersample_task_tracking_verified;
        stats.intersample_direction_scope='stored accepted segments only; gaps not bridged';
        trajectory.intersample_direction_scope=stats.intersample_direction_scope;
        stats.verification_method=[stats.verification_method,'; adaptive whole-cell direction envelope'];
    end
    % Whole-input feasibility additionally requires every input row covered.
    if stats.n_ok>0 && ~stats.verification.passed && ~coverage_only
        error('coneik_inverse_trajectory:InternalVerification','Accepted output failed its final coefficient verification.');
    end
end

function [segment,entry]=solve_segment(t,gamma,cache,rows,opts,total_rows)
    segment=struct();n=numel(t);
    entry=struct('rows',rows,'status','no_feasible_candidate_path','failure_reason','',...
        'candidate_evaluations',0,'necessary_speed_pruned',0,'beam_pruned',0,...
        'max_layer_candidates',0,'max_layer_retained',0,'candidate_limit_hit',false,...
        'center_family_sampled',false,'periodic_lifts_limited',false,...
        'position_enumeration_complete',true,'retained_full_paths',0,...
        'verified_full_paths',0,'passing_full_paths',0,'unchecked_full_paths',0,...
        'candidate_diagnostics_complete',false,'selected_cost',NaN,...
        'selected_path_index',NaN,'selected_branches',[],'candidate_costs',[],...
        'candidate_status',{{}},'candidate_verification',{{}},'verification',struct(),...
        'search_complete',false,'global_feasibility_checked',false,...
        'intersample_task_tracking_verified',false,...
        'scope','finite representatives and displacement-ranked beam; no global feasibility claim');
    boundary_enabled=any(isfield(opts,{'initial_q','initial_velocity','final_velocity','joint_excursion_bounds'}));
    effective_bounds=opts.joint_bounds;bounds_ok=true;
    if isfield(opts,'joint_excursion_bounds')
        window=opts.initial_q+opts.joint_excursion_bounds;
        bounds_ok=all(isfinite(window(:)))&&...
            ~any(opts.joint_excursion_bounds(:)~=0 & window(:)==repmat(opts.initial_q,2,1));
        if bounds_ok,effective_bounds=[max(effective_bounds(:,1),window(:,1)),min(effective_bounds(:,2),window(:,2))];end
    end
    verify_opts=opts;
    if boundary_enabled
        verify_opts.joint_bounds=effective_bounds;
        if isfield(verify_opts,'joint_excursion_bounds'),verify_opts=rmfield(verify_opts,'joint_excursion_bounds');end
        if rows(1)~=1
            for name={'initial_q','initial_velocity'}
                if isfield(verify_opts,name{1}),verify_opts=rmfield(verify_opts,name{1});end
            end
        end
        if rows(end)~=total_rows&&isfield(verify_opts,'final_velocity'),verify_opts=rmfield(verify_opts,'final_velocity');end
    end
    if n<2
        entry.status='insufficient_segment_samples';
        entry.failure_reason='A singleton does not determine continuous motion.';return;
    end
    if ~bounds_ok
        entry.status='search_unresolved';entry.failure_reason='excursion_absolute_window_unresolved';return;
    end
    hard_start=isfield(verify_opts,'initial_q');start_branch=NaN;
    if hard_start
        [start_branch,start_status,start_reason]=check_initial_q(opts.initial_q,gamma(1,:)',cache{1},effective_bounds,opts);
        if ~strcmp(start_status,'pass'),entry.status=start_status;entry.failure_reason=start_reason;return;end
    end
    if any(effective_bounds(:,1)>effective_bounds(:,2))
        entry.failure_reason='empty_effective_joint_bounds';return;
    end
    qlayer=cell(n,1);player=cell(n,1);blayer=cell(n,1);cost=[];
    resolved_model=[];if isfield(opts,'model'),resolved_model=opts.model;end
    position_opts=struct('joint_bounds',effective_bounds,'max_candidates',opts.max_candidates,...
        'winding_radius',opts.winding_radius,'center_q1',opts.center_q1);
    if isfield(opts,'model'),position_opts.model=opts.model;end
    for i=1:n
        if i==1,parents=1;else,parents=numel(cost);end
        qchunks=cell(1,parents);parentchunks=cell(1,parents);
        branchchunks=cell(1,parents);costchunks=cell(1,parents);
        % Only a complete finite non-center set is reference-independent.
        % Unbounded lifts and center seeds still depend on each parent.
        reuse_candidates=all(isfinite(effective_bounds(:)))&&...
            ~strcmp(cache{i}.position_info.status,'center')&&~(i==1&&hard_start);
        for p=1:parents
            if i==1
                reference=opts.previous_q;
            else
                reference=qlayer{i-1}(:,p);
            end
            position_opts.previous_q=reference;
            try
                if i==1&&hard_start
                    Q=opts.initial_q;
                    meta=struct('branch',start_branch,'base_candidates',cache{i}.base_candidates,...
                        'center_family_sampled',false,'periodic_lifts_limited',false,'position_enumeration_complete',true);
                elseif p==1||~reuse_candidates
                    [Q,meta]=coneik_position_candidates_resolved(gamma(i,:)',position_opts,cache{i},resolved_model,false);
                end
            catch exception
                entry.failure_reason=exception.identifier;
                if strcmp(exception.identifier,'coneik_position_candidates:CandidateLimitExceeded')
                    entry.status='search_budget_exceeded';entry.candidate_limit_hit=true;
                else
                    entry.status='search_unresolved';
                end
                return;
            end
            entry.candidate_evaluations=entry.candidate_evaluations+size(Q,2);
            entry.center_family_sampled=entry.center_family_sampled||meta.center_family_sampled;
            entry.periodic_lifts_limited=entry.periodic_lifts_limited||meta.periodic_lifts_limited;
            entry.position_enumeration_complete=entry.position_enumeration_complete&&meta.position_enumeration_complete;
            if i==1
                % Without an explicit reference, a branch's base representative
                % has no arbitrary zero-position penalty. Prefer fewer turns.
                if isempty(reference)
                    anchor=zeros(size(Q));
                    for j=1:size(Q,2)
                        b=meta.branch(j);if b==0,b=1;end
                        anchor(:,j)=meta.base_candidates(:,b);
                    end
                    if meta.center_family_sampled,anchor(1,:)=Q(1,:);end
                    score=sum((Q-anchor).^2,1);
                else
                    score=sum((Q-reference).^2,1);
                end
                parentids=zeros(1,size(Q,2));keep=true(1,size(Q,2));
            else
                delta=Q-reference;dt=t(i)-t(i-1);
                speedtol=1e-8;
                allowed=opts.velocity_limits*dt+speedtol*dt;
                keep=all(abs(delta)<=allowed,1);
                entry.necessary_speed_pruned=entry.necessary_speed_pruned+nnz(~keep);
                score=cost(p)+sum(delta.^2,1);
                parentids=p*ones(1,size(Q,2));
            end
            if any(~isfinite(score(keep)))
                entry.status='search_unresolved';entry.failure_reason='nonfinite_path_cost';return;
            end
            qchunks{p}=Q(:,keep);parentchunks{p}=parentids(keep);
            branchchunks{p}=meta.branch(keep);costchunks{p}=score(keep);
        end
        % Preserve parent order without repeatedly copying the whole layer.
        nextq=[qchunks{:}];nextparent=[parentchunks{:}];
        nextbranch=[branchchunks{:}];nextcost=[costchunks{:}];
        count=numel(nextcost);entry.max_layer_candidates=max(entry.max_layer_candidates,count);
        if count==0
            entry.failure_reason='No candidate survives position bounds and necessary speed-edge checks.';return;
        end
        % Ties use branch, physical angle, then parent, giving reproducible paths.
        [~,order]=sortrows([nextcost(:),nextbranch(:),nextq',nextparent(:)],[1 2 3 4 5]);
        retained=order(1:min(opts.beam_width,count));
        entry.beam_pruned=entry.beam_pruned+count-numel(retained);
        qlayer{i}=nextq(:,retained);player{i}=nextparent(retained);blayer{i}=nextbranch(retained);
        cost=nextcost(retained);entry.max_layer_retained=max(entry.max_layer_retained,numel(retained));
    end
    entry.retained_full_paths=numel(cost);entry.candidate_costs=cost;
    entry.candidate_status=repmat({'not_checked'},1,numel(cost));
    entry.candidate_verification=cell(1,numel(cost));
    entry.unchecked_full_paths=numel(cost);
    for p=1:numel(cost)
        q=zeros(n,2);branches=zeros(n,1);index=p;
        for i=n:-1:1
            q(i,:)=qlayer{i}(:,index)';branches(i)=blayer{i}(index);index=player{i}(index);
        end
        try
            pp={};interpolation='not-a-knot spline';
            if isfield(verify_opts,'initial_velocity')||isfield(verify_opts,'final_velocity')
                v0=[];vf=[];
                if isfield(verify_opts,'initial_velocity'),v0=verify_opts.initial_velocity;end
                if isfield(verify_opts,'final_velocity'),vf=verify_opts.final_velocity;end
                for axis=1:2
                    left=[];right=[];if ~isempty(v0),left=v0(axis);end;if ~isempty(vf),right=vf(axis);end
                    pp{axis}=coneik_boundary_spline(t,q(:,axis),left,right); %#ok<AGROW>
                end
                interpolation='C2 cubic with requested endpoint velocity; free side not-a-knot (two points quadratic)';
            else
                pp={spline(t,q(:,1)),spline(t,q(:,2))};
            end
            candidate=struct('rows',rows,'t',t,'gamma',gamma,...
                'pp',{pp},'branches',branches,'q_samples',q,...
                'path_cost',cost(p),'interpolation',interpolation);
            if isfield(opts,'model'),candidate.model_key=opts.model.key;end
            verification=coneik_verify_trajectory(struct('segments',candidate),verify_opts);
        catch exception
            verification=coneik_verify_trajectory(struct('segments',struct([])),verify_opts);
            verification.reason=['spline_construction_unresolved: ',exception.identifier];
        end
        entry.verified_full_paths=entry.verified_full_paths+1;
        entry.unchecked_full_paths=numel(cost)-entry.verified_full_paths;
        entry.candidate_status{p}=verification.status;
        entry.candidate_verification{p}=verification;
        if verification.passed
            entry.passing_full_paths=entry.passing_full_paths+1;
            if isempty(fieldnames(segment))
                segment=candidate;entry.status='pass';entry.failure_reason='';
                entry.selected_cost=cost(p);entry.selected_path_index=p;
                entry.selected_branches=branches;entry.verification=verification;
                if isfield(opts,'intersample_direction_tolerance_deg')
                    entry.intersample_task_tracking_verified=verification.intersample_task_tracking_verified;
                end
            end
            if ~opts.full_candidate_diagnostics,break;end
        end
    end
    entry.candidate_diagnostics_complete=entry.unchecked_full_paths==0;
    if ~strcmp(entry.status,'pass')
        entry.failure_reason='No retained full path passed continuous polynomial joint bounds and sample FK/C2 checks.';
        if isfield(opts,'intersample_direction_tolerance_deg')
            entry.failure_reason='No retained full path passed joint bounds, sample FK/C2 and intersample direction checks.';
        end
        if any(strcmp(entry.candidate_status,'unresolved'))
            entry.status='search_unresolved';
            entry.failure_reason='At least one retained full path could not be numerically verified; none passed.';
        end
    end
end

function opts=validate_options(options)
    allowed={'mode','joint_bounds','velocity_limits','acceleration_limits','previous_q',...
        'beam_width','max_candidates','winding_radius','center_q1','model','position_backend',...
        'intersample_direction_tolerance_deg','intersample_direction_max_refinements',...
        'initial_q','initial_velocity','final_velocity','joint_excursion_bounds','full_candidate_diagnostics'};
    if isstruct(options)&&isscalar(options)
        if any(isfield(options,{'trajectory_selection_policy','trajectory_selection_weights'}))
            error('coneik_inverse_trajectory:RemovedFeature',...
                'Trajectory smoothness ranking is temporarily removed. Omit trajectory_selection_policy and trajectory_selection_weights; the first verified candidate is used.');
        end
    end
    if ~(isstruct(options)&&isscalar(options)) || any(~ismember(fieldnames(options),allowed))
        error('coneik_inverse_trajectory:InvalidOptions','Unknown option or non-scalar struct.');
    end
    opts=struct('mode','continuous','joint_bounds',[-Inf Inf;-Inf Inf],...
        'velocity_limits',[Inf;Inf],'acceleration_limits',[Inf;Inf],...
        'previous_q',[],'beam_width',8,'max_candidates',256,'winding_radius',1,'center_q1',[],...
        'full_candidate_diagnostics',false);
    names=fieldnames(options);for k=1:numel(names),opts.(names{k})=options.(names{k});end
    if ~(islogical(opts.full_candidate_diagnostics)&&isscalar(opts.full_candidate_diagnostics))
        error('coneik_inverse_trajectory:InvalidOptions','full_candidate_diagnostics must be a logical scalar.');
    end
    for name={'initial_q','initial_velocity','final_velocity'}
        if isfield(opts,name{1})
            value=opts.(name{1});
            if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&all(isfinite(value(:))))
                error('coneik_inverse_trajectory:InvalidOptions','Boundary q/v must be finite real 2-vectors.');
            end
            opts.(name{1})=double(value(:));
        end
    end
    if isfield(opts,'joint_excursion_bounds')
        b=opts.joint_excursion_bounds;
        if ~isfield(opts,'initial_q')||~(isnumeric(b)&&isreal(b)&&isequal(size(b),[2 2])&&all(isfinite(b(:)))&&all(b(:,1)<=0)&&all(b(:,2)>=0))
            error('coneik_inverse_trajectory:InvalidOptions','Excursion requires initial_q and finite signed [lower<=0 upper>=0] rows.');
        end
        opts.joint_excursion_bounds=double(b);
    end
    if isfield(opts,'intersample_direction_tolerance_deg')
        v=opts.intersample_direction_tolerance_deg;
        if ~(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v)&&v>0&&v<180)
            error('coneik_inverse_trajectory:InvalidOptions','Direction tolerance must satisfy 0<degrees<180.');
        end
        opts.intersample_direction_tolerance_deg=double(v);
        chord=2*sind(double(v)/2);
        if ~(isfinite(chord)&&chord>0)
            error('coneik_inverse_trajectory:InvalidOptions','Direction tolerance is not representable.');
        end
        if ~isfield(opts,'intersample_direction_max_refinements'),opts.intersample_direction_max_refinements=1024;end
        v=opts.intersample_direction_max_refinements;
        if ~(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v)&&v>=0&&fix(v)==v&&v<=flintmax/2)
            error('coneik_inverse_trajectory:InvalidOptions','Direction refinement budget must be a nonnegative integer.');
        end
        opts.intersample_direction_max_refinements=double(v);
    elseif isfield(opts,'intersample_direction_max_refinements')
        error('coneik_inverse_trajectory:InvalidOptions','A direction refinement budget requires an explicit tolerance.');
    end
    if isfield(opts,'model'),opts.model=coneik_resolve_model(opts.model);end
    if isfield(opts,'position_backend')
        backend=opts.position_backend;
        if isstring(backend)&&isscalar(backend)&&~ismissing(backend),backend=char(backend);end
        if ~(ischar(backend)&&isrow(backend)&&any(strcmp(backend,{'matlab','cpp','auto'})))
            error('coneik_inverse_trajectory:InvalidOptions','position_backend must be matlab, cpp or auto.');
        end
        opts.position_backend=backend;
    end
    if ~((ischar(opts.mode)&&isrow(opts.mode))||(isstring(opts.mode)&&isscalar(opts.mode)&&~ismissing(opts.mode))) || ~strcmp(opts.mode,'continuous')
        error('coneik_inverse_trajectory:InvalidOptions','mode must be continuous.');
    end
    b=opts.joint_bounds;
    if ~(isnumeric(b)&&isreal(b)&&isequal(size(b),[2 2])&&~any(isnan(b(:)))&&all(b(:,1)<=b(:,2))&&~any(b(:,1)==Inf)&&~any(b(:,2)==-Inf))
        error('coneik_inverse_trajectory:InvalidOptions','joint_bounds needs valid [lower upper] rows.');
    end
    opts.joint_bounds=double(b);
    for name={'velocity_limits','acceleration_limits'}
        value=opts.(name{1});
        if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&all(value(:)>=0))
            error('coneik_inverse_trajectory:InvalidOptions','Limits must be nonnegative 2-vectors.');
        end
        opts.(name{1})=double(value(:));
    end
    if ~isempty(opts.previous_q)
        if ~(isnumeric(opts.previous_q)&&isreal(opts.previous_q)&&isvector(opts.previous_q)&&numel(opts.previous_q)==2&&all(isfinite(opts.previous_q(:))))
            error('coneik_inverse_trajectory:InvalidOptions','previous_q must be a finite 2-vector.');
        end
        opts.previous_q=double(opts.previous_q(:));
    end
    if ~(isnumeric(opts.center_q1)&&isreal(opts.center_q1)&&(isvector(opts.center_q1)||isempty(opts.center_q1))&&all(isfinite(opts.center_q1(:))))
        error('coneik_inverse_trajectory:InvalidOptions','center_q1 must be finite seeds.');
    end
    for name={'beam_width','max_candidates','winding_radius'}
        value=opts.(name{1});minimum=1;if strcmp(name{1},'winding_radius'),minimum=0;end
        if ~(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=minimum&&fix(value)==value&&value<=flintmax/2)
            error('coneik_inverse_trajectory:InvalidOptions','Invalid integer search budget.');
        end
    end
end

function [branch,status,reason]=check_initial_q(q,gamma,cache,bounds,opts)
    branch=NaN;status='boundary_incompatible';reason='initial_q_outside_hard_bounds';
    if any(q<bounds(:,1)-1e-9|q>bounds(:,2)+1e-9),return;end
    turns=q/360;
    if any(~isfinite(turns)|abs(turns)>flintmax/2)
        status='search_unresolved';reason='initial_q_periodic_resolution_unresolved';return;
    end
    try
        model=[];if isfield(opts,'model'),model=opts.model;end
        r=coneik_kinematics(q,[0;0],model);
        target=[cosd(gamma(1))*cosd(gamma(2));cosd(gamma(1))*sind(gamma(2));-sind(gamma(1))];
        residual=norm(r-target);
        if ~isfinite(residual),status='search_unresolved';reason='initial_q_fk_unresolved';return;end
    catch
        status='search_unresolved';reason='initial_q_fk_unresolved';return;
    end
    reason='initial_q_target_mismatch';if residual>1e-10,return;end
    kind=cache.position_info.status;
    if strcmp(kind,'center'),branch=0;
    elseif any(strcmp(kind,{'inner_boundary','outer_boundary'})),branch=1;
    else
        for j=1:2
            base=cache.base_candidates(:,j);turns=round((q-base)/360);
            if any(~isfinite(turns)|abs(turns)>flintmax/2)
                status='search_unresolved';reason='initial_q_periodic_resolution_unresolved';return;
            end
            reconstructed=base+360*turns;
            if any(~isfinite(reconstructed))
                status='search_unresolved';reason='initial_q_periodic_resolution_unresolved';return;
            end
            if all(abs(reconstructed-q)<=1e-9),branch=j;break;end
        end
    end
    if isnan(branch),reason='initial_q_not_on_analytic_branch';return;end
    status='pass';reason='';
end
