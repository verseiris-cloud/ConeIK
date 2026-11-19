function report = coneik_verify_trajectory(trajectory, options, result)
%CONEIK_VERIFY_TRAJECTORY Independently check stored cubic joint trajectories.
% report=coneik_verify_trajectory(trajectory,options,result_optional).
% OPTIONS: joint_bounds is 2x2 [lower upper] in nominal degrees;
% velocity_limits / acceleration_limits are nonnegative 2-vectors.
% Valid numeric limits are converted to double before bound arithmetic.
% model is optional: a stored trajectory model is used when omitted.
% An explicit model must match stored model metadata, including segment keys.
% Endpoints and stationary roots cover the entire polynomial domains.
% Float64 checks use stated engineering tolerances, not interval arithmetic.
% Sample FK and C2 joins are always checked. Presence of
% intersample_direction_tolerance_deg (0<degrees<180) additionally checks
% shortest-arc SLERP between sample directions, linear in physical time.
% intersample_direction_max_refinements defaults to 1024 bisections/segment.
% This optional float64 bound is not an interval-arithmetic certificate.
% Full documented trajectory option structs are accepted for convenience;
% full_candidate_diagnostics is accepted as a search-diagnostic no-op;
% it cannot disable any verification. Unknown names are errors. Extreme polynomial scales may be unresolved.
% Optional initial_q / initial_velocity / final_velocity apply at the original
% source_t endpoints (or stored-domain endpoints for external trajectories).
% joint_excursion_bounds is a signed physical window from initial_q, checked
% over every stored cubic. Endpoint acceleration is diagnostic, not constrained.
    if nargin<2 || isempty(options), options=struct(); end
    if ~(isstruct(options)&&isscalar(options))
        error('coneik_verify_trajectory:InvalidOptions','options must be a scalar struct.');
    end
    allowed={'mode','joint_bounds','velocity_limits','acceleration_limits',...
        'previous_q','beam_width','max_candidates','winding_radius','center_q1','model','position_backend',...
        'intersample_direction_tolerance_deg','intersample_direction_max_refinements',...
        'initial_q','initial_velocity','final_velocity','joint_excursion_bounds','full_candidate_diagnostics'};
    if any(isfield(options,{'trajectory_selection_policy','trajectory_selection_weights'}))
        error('coneik_verify_trajectory:RemovedFeature',...
            'Trajectory smoothness ranking is temporarily removed. Omit trajectory_selection_policy and trajectory_selection_weights; the first verified candidate is used.');
    end
    if any(~ismember(fieldnames(options),allowed))
        error('coneik_verify_trajectory:InvalidOptions','Unknown verifier option field.');
    end
    if isfield(options,'full_candidate_diagnostics') && ...
            ~(islogical(options.full_candidate_diagnostics)&&isscalar(options.full_candidate_diagnostics))
        error('coneik_verify_trajectory:InvalidOptions','full_candidate_diagnostics must be a logical scalar.');
    end
    bounds=get_option(options,'joint_bounds',[-Inf Inf;-Inf Inf]);
    velocity=get_option(options,'velocity_limits',[Inf;Inf]);
    acceleration=get_option(options,'acceleration_limits',[Inf;Inf]);
    validate_limits(bounds,velocity,acceleration);
    tracking=tracking_options(options);
    boundary=boundary_options(options);
    bounds=double(bounds);
    velocity=double(velocity(:)'); acceleration=double(acceleration(:)');
    tol=struct('position',[1e-9 1e-9],...
        'velocity',[1e-8 1e-8],'acceleration',[1e-7 1e-7],...
        'sample_direction_chord',1e-10,'sample_joint',[1e-9 1e-9]);
    report=struct('status','unresolved','passed',false,'reason','',...
        'method','float64 polynomial endpoints and stationary roots',...
        'certified',false,'global_feasibility_checked',false,...
        'continuous_joint_bounds_checked',false,...
        'intersample_task_tracking_verified',false,'initial_velocity','free',...
        'tolerance',tol,'joint_bounds',bounds,'velocity_limits',velocity,...
        'acceleration_limits',acceleration,'segment_count',0,'segments',struct([]),...
        'sample_direction_residual_max',NaN,'result_consistency_max',NaN);
    if boundary.enabled
        report.boundary_constraints=empty_boundary(boundary);
        if boundary.has_initial_velocity,report.initial_velocity='hard_requested';end
    end
    if tracking.enabled
        report.intersample_direction=empty_tracking(tracking);
        report.intersample_direction_scope='stored accepted segments only; gaps not bridged';
    end
    if ~(isstruct(trajectory)&&isscalar(trajectory)&&isfield(trajectory,'segments'))
        report.reason='invalid_trajectory'; return;
    end
    model=verification_model(trajectory,options);
    if ~isempty(model),report.model_key=model.key;end
    if isempty(trajectory.segments)
        report.reason='no_accepted_segments';
        if boundary.enabled
            report.boundary_constraints=verify_boundaries(trajectory,report.segments,boundary,tol,bounds);
            report.reason=report.boundary_constraints.reason;
        end
        return;
    end
    report.segment_count=numel(trajectory.segments);
    statuses=cell(1,report.segment_count);
    for s=1:report.segment_count
        try
            part=verify_segment(trajectory.segments(s),bounds,velocity,acceleration,tol,model);
        catch exception
            part=empty_part(); part.reason=['invalid_polynomial: ',exception.identifier];
        end
        if tracking.enabled
            part.intersample_direction=empty_tracking(tracking);
            if part.passed
                part.intersample_direction=verify_direction(trajectory.segments(s),part.extrema,tracking,model);
                part.status=part.intersample_direction.status;
                part.passed=strcmp(part.status,'pass');
                if ~part.passed,part.reason=['intersample_direction: ',part.intersample_direction.reason];end
            end
        end
        if s==1,report.segments=part;else,report.segments(s)=part;end
        statuses{s}=part.status;
    end
    report.sample_direction_residual_max=max([report.segments.sample_direction_residual_max]);
    report.continuous_joint_bounds_checked=all([report.segments.extrema_resolved]);
    if any(strcmp(statuses,'unresolved'))
        report.reason='one_or_more_segments_unresolved';
    elseif any(strcmp(statuses,'violation'))
        report.status='violation'; report.reason='one_or_more_checks_failed';
    else
        report.status='pass'; report.passed=true;
    end
    domains=zeros(report.segment_count,2);
    for s=1:report.segment_count
        if isfield(trajectory.segments(s),'t') && numel(trajectory.segments(s).t)>=2
            ts=trajectory.segments(s).t;domains(s,:)=[ts(1) ts(end)];
        else,domains(s,:)=[NaN NaN];end
    end
    domains=sortrows(domains,1);
    if any(domains(2:end,1)<=domains(1:end-1,2))
        report.status='violation';report.passed=false;report.reason='overlapping_segment_domains';
    end
    if nargin>=3
        try
            if ~(isnumeric(result)&&isreal(result)&&size(result,2)==7&&all(isfinite(result(:,1))))
                error('coneik_verify_trajectory:InvalidResult','Invalid result matrix.');
            end
            expected=coneik_trajectory_eval(trajectory,result(:,1));
            if ~isequal(isnan(expected(:,2:7)),isnan(result(:,2:7))) || any(isinf(result(:)))
                report.status='violation';report.passed=false;report.reason='result_domain_or_finiteness_mismatch';
            else
                delta=abs(expected(:,2:7)-result(:,2:7));
                valid=isfinite(delta); vals=delta(valid);
                if isempty(vals),report.result_consistency_max=0;else,report.result_consistency_max=max(vals);end
                threshold=[tol.position tol.velocity tol.acceleration];
                if any(any(delta>threshold))
                    report.status='violation';report.passed=false;report.reason='result_not_from_stored_coefficients';
                end
            end
        catch exception
            report.status='unresolved';report.passed=false;report.reason=['result_check_failed: ',exception.identifier];
        end
    end
    if boundary.enabled
        checks=verify_boundaries(trajectory,report.segments,boundary,tol,bounds);
        checks.coverage_only_failure=checks.coverage_only_failure&&report.passed;
        if report.passed&&~checks.verified
            report.status=checks.status;report.passed=false;report.reason=checks.reason;
        end
        report.boundary_constraints=checks;
    end
    if tracking.enabled
        report.intersample_direction=aggregate_tracking(report.segments,tracking);
        if boundary.enabled&&report.boundary_constraints.coverage_only_failure
            % Missing source endpoints do not erase verified stored-segment tracking.
            report.intersample_task_tracking_verified=strcmp(report.intersample_direction.status,'pass');
        else
            report.intersample_task_tracking_verified=report.passed&&strcmp(report.intersample_direction.status,'pass');
        end
    end
end

function options=boundary_options(input)
    names={'initial_q','initial_velocity','final_velocity','joint_excursion_bounds'};
    options=struct('enabled',any(isfield(input,names)),...
        'has_initial_q',false,'has_initial_velocity',false,'has_final_velocity',false,'has_excursion',false,...
        'initial_q',[],'initial_velocity',[],'final_velocity',[],'joint_excursion_bounds',[]);
    for k=1:3
        name=names{k};
        if isfield(input,name)
            v=input.(name);
            if ~(isnumeric(v)&&isreal(v)&&isvector(v)&&numel(v)==2&&all(isfinite(v(:))))
                error('coneik_verify_trajectory:InvalidOptions','Boundary q/v must be finite real 2-vectors.');
            end
            options.(name)=double(v(:));options.(['has_',name])=true;
        end
    end
    if isfield(input,'joint_excursion_bounds')
        v=input.joint_excursion_bounds;
        if ~options.has_initial_q||~(isnumeric(v)&&isreal(v)&&isequal(size(v),[2 2])&&all(isfinite(v(:)))&&all(v(:,1)<=0)&&all(v(:,2)>=0))
            error('coneik_verify_trajectory:InvalidOptions','Excursion requires initial_q and finite signed [lower<=0 upper>=0] rows.');
        end
        options.has_excursion=true;options.joint_excursion_bounds=double(v);
    end
end

function out=empty_boundary(b)
    out=struct('requested_initial_q',b.has_initial_q,'requested_initial_velocity',b.has_initial_velocity,...
        'requested_final_velocity',b.has_final_velocity,'excursion_requested',b.has_excursion,...
        'source_start_time',NaN,'source_end_time',NaN,'boundary_time_source','unavailable',...
        'initial_endpoint_covered',false,'final_endpoint_covered',false,...
        'initial_q_residual',[NaN NaN],'initial_velocity_residual',[NaN NaN],...
        'final_velocity_residual',[NaN NaN],'observed_initial_acceleration',[NaN NaN],...
        'observed_final_acceleration',[NaN NaN],'excursion_relative_min',[NaN NaN],...
        'excursion_relative_max',[NaN NaN],'excursion_margins',nan(2,2),...
        'effective_joint_bounds',nan(2,2),'segments',struct([]),...
        'acceleration_boundary_constrained',false,'gap_bridged',false,...
        'verified',false,'status','unresolved','reason','boundary_domain_unavailable','coverage_only_failure',false);
end

function out=verify_boundaries(trajectory,parts,b,tol,bounds)
    out=empty_boundary(b);out.effective_joint_bounds=bounds;
    if isfield(trajectory,'source_t')
        times=trajectory.source_t;
        if ~(isnumeric(times)&&isreal(times)&&isvector(times)&&~isempty(times)&&all(isfinite(times(:)))&&all(isfinite(diff(double(times(:)))))&&all(diff(double(times(:)))>0))
            out.reason='boundary_source_domain_unresolved';return;
        end
        first=double(times(1));last=double(times(end));out.boundary_time_source='source_t';
    else
        domains=zeros(numel(trajectory.segments),2);
        for k=1:numel(trajectory.segments)
            seg=trajectory.segments(k);
            if ~isfield(seg,'t')||~isnumeric(seg.t)||~isreal(seg.t)||numel(seg.t)<2||any(~isfinite(seg.t(:)))
                out.reason='boundary_stored_domain_unresolved';return;
            end
            domains(k,:)=[seg.t(1),seg.t(end)];
        end
        if isempty(domains),out.reason='boundary_domain_unavailable';return;end
        first=min(domains(:,1));last=max(domains(:,2));out.boundary_time_source='stored_segment_domain';
    end
    out.source_start_time=first;out.source_end_time=last;
    if b.has_excursion
        window=b.initial_q+b.joint_excursion_bounds;
        if any(~isfinite(window(:)))||any(b.joint_excursion_bounds(:)~=0&window(:)==repmat(b.initial_q,2,1))
            out.reason='excursion_absolute_window_unresolved';return;
        end
        out.effective_joint_bounds=[max(bounds(:,1),window(:,1)),min(bounds(:,2),window(:,2))];
    end
    relative_min=[Inf Inf];relative_max=[-Inf -Inf];
    for k=1:numel(trajectory.segments)
        seg=trajectory.segments(k);
        item=struct('status','unresolved','passed',false,'reason','boundary_segment_unresolved',...
            'initial_applies',false,'final_applies',false,'endpoint_q',nan(2,2),...
            'endpoint_velocity',nan(2,2),'endpoint_acceleration',nan(2,2),...
            'initial_q_residual',[NaN NaN],'initial_velocity_residual',[NaN NaN],...
            'final_velocity_residual',[NaN NaN],'excursion_relative_min',[NaN NaN],...
            'excursion_relative_max',[NaN NaN],'excursion_margins',nan(2,2));
        try
            if ~isfield(seg,'t')||numel(seg.t)<2||seg.t(1)<first||seg.t(end)>last
                item.status='violation';item.reason='stored_segment_outside_source_domain';
            else
                item.initial_applies=seg.t(1)==first;item.final_applies=seg.t(end)==last;
                out.initial_endpoint_covered=out.initial_endpoint_covered||item.initial_applies;
                out.final_endpoint_covered=out.final_endpoint_covered||item.final_applies;
                for axis=1:2
                    [breaks,c,~,order,dim]=unmkpp(seg.pp{axis});
                    if dim~=1||order>4||order<1||any(~isfinite(c(:)))||breaks(1)~=seg.t(1)||breaks(end)~=seg.t(end)
                        error('coneik_verify_trajectory:NumericalUnresolved','Invalid boundary polynomial.');
                    end
                    c=[zeros(size(c,1),4-order),c];h=breaks(end)-breaks(end-1);
                    vc=[3*c(:,1),2*c(:,2),c(:,3)];ac=[6*c(:,1),2*c(:,2)];
                    item.endpoint_q(:,axis)=[c(1,4);polyval(c(end,:),h)];
                    item.endpoint_velocity(:,axis)=[vc(1,end);polyval(vc(end,:),h)];
                    item.endpoint_acceleration(:,axis)=[ac(1,end);polyval(ac(end,:),h)];
                end
                if any(~isfinite([item.endpoint_q(:);item.endpoint_velocity(:);item.endpoint_acceleration(:)]))
                    error('coneik_verify_trajectory:NumericalUnresolved','Nonfinite one-sided boundary derivative.');
                end
                item.status='pass';item.passed=true;item.reason='';
                if item.initial_applies
                    out.observed_initial_acceleration=item.endpoint_acceleration(1,:);
                    if b.has_initial_q
                        item.initial_q_residual=abs(boundary_difference(item.endpoint_q(1,:),b.initial_q'));
                        out.initial_q_residual=item.initial_q_residual;
                        if any(item.initial_q_residual>tol.position),item=boundary_violation(item,'initial_q_mismatch');end
                    end
                    if b.has_initial_velocity
                        item.initial_velocity_residual=abs(boundary_difference(item.endpoint_velocity(1,:),b.initial_velocity'));
                        out.initial_velocity_residual=item.initial_velocity_residual;
                        if any(item.initial_velocity_residual>tol.velocity),item=boundary_violation(item,'initial_velocity_mismatch');end
                    end
                end
                if item.final_applies
                    out.observed_final_acceleration=item.endpoint_acceleration(2,:);
                    if b.has_final_velocity
                        item.final_velocity_residual=abs(boundary_difference(item.endpoint_velocity(2,:),b.final_velocity'));
                        out.final_velocity_residual=item.final_velocity_residual;
                        if any(item.final_velocity_residual>tol.velocity),item=boundary_violation(item,'final_velocity_mismatch');end
                    end
                end
                if b.has_excursion
                    if k>numel(parts)||~parts(k).extrema_resolved
                        error('coneik_verify_trajectory:NumericalUnresolved','Excursion extrema unresolved.');
                    end
                    lo=boundary_difference(parts(k).extrema.q_min,b.initial_q');
                    hi=boundary_difference(parts(k).extrema.q_max,b.initial_q');
                    item.excursion_relative_min=lo;item.excursion_relative_max=hi;
                    margins=[boundary_difference(lo',b.joint_excursion_bounds(:,1)),boundary_difference(b.joint_excursion_bounds(:,2),hi')];
                    item.excursion_margins=margins;
                    relative_min=min(relative_min,lo);relative_max=max(relative_max,hi);
                    if any(any(margins< -tol.position(:))),item=boundary_violation(item,'joint_excursion_violation');end
                end
            end
        catch exception
            item.status='unresolved';item.passed=false;item.reason=['boundary_arithmetic_unresolved: ',exception.identifier];
        end
        if k==1,out.segments=item;else,out.segments(k)=item;end
    end
    if b.has_excursion&&all(isfinite([relative_min relative_max]))
        out.excursion_relative_min=relative_min;out.excursion_relative_max=relative_max;
        out.excursion_margins=[relative_min'-b.joint_excursion_bounds(:,1),b.joint_excursion_bounds(:,2)-relative_max'];
    end
    if ~isempty(out.segments)&&any(strcmp({out.segments.status},'violation'))
        out.status='violation';out.reason='one_or_more_boundary_checks_failed';return;
    end
    if ~isempty(out.segments)&&any(strcmp({out.segments.status},'unresolved'))
        out.reason='one_or_more_boundary_checks_unresolved';return;
    end
    missing_initial=(b.has_initial_q||b.has_initial_velocity)&&~out.initial_endpoint_covered;
    missing_final=b.has_final_velocity&&~out.final_endpoint_covered;
    if missing_initial||missing_final
        if missing_initial&&missing_final,out.reason='requested_initial_and_final_boundary_not_covered';
        elseif missing_initial,out.reason='requested_initial_boundary_not_covered';
        else,out.reason='requested_final_boundary_not_covered';end
        out.coverage_only_failure=~isempty(out.segments)&&all([out.segments.passed]);return;
    end
    if isempty(out.segments),out.reason='no_accepted_segments';return;end
    out.status='pass';out.verified=true;out.reason='';
end
function value=boundary_difference(left,right)
    value=left-right;
    if any(~isfinite(value(:)))||any(left(:)~=right(:)&(value(:)==0|abs(value(:))<realmin))
        error('coneik_verify_trajectory:NumericalUnresolved','Boundary subtraction is not representable.');
    end
end
function item=boundary_violation(item,reason)
    item.status='violation';item.passed=false;item.reason=reason;
end

function part=verify_segment(seg,bounds,vlimit,alimit,tol,model)
    part=empty_part();
    if ~(isfield(seg,'t')&&isfield(seg,'gamma')&&isfield(seg,'pp')&&numel(seg.pp)==2)
        part.reason='missing_segment_fields';return;
    end
    t=seg.t(:); n=numel(t);
    if n<2 || ~isnumeric(t) || ~isreal(t) || any(~isfinite(t)) || any(~isfinite(diff(t))) || any(diff(t)<=0) ||...
            ~isnumeric(seg.gamma) || ~isreal(seg.gamma) || ~isequal(size(seg.gamma),[n 2]) || any(~isfinite(seg.gamma(:)))
        part.reason='invalid_segment_samples';return;
    end
    extrema=struct('q_min',zeros(1,2),'q_max',zeros(1,2),'velocity_abs_max',zeros(1,2),'acceleration_abs_max',zeros(1,2));
    when=extrema; q=zeros(n,2); joinmax=zeros(1,3);
    endpoint_observed=zeros(2,3);endpoint_discrepancy=zeros(2,3);endpoint_tolerance=zeros(2,3);
    for axis=1:2
        [breaks,c,pieces,order,dim]=unmkpp(seg.pp{axis});
        if dim~=1 || order>4 || order<1 || ~isreal(c) || ~isreal(breaks) || any(~isfinite(c(:))) || any(~isfinite(breaks)) || any(~isfinite(diff(breaks))) || any(diff(breaks)<=0)
            part.reason='nonfinite_or_invalid_polynomial';return;
        end
        if breaks(1)~=t(1) || breaks(end)~=t(end)
            part.reason='polynomial_domain_mismatch';return;
        end
        c=[zeros(pieces,4-order),c]; h=diff(breaks(:));
        [P,VP,AP,original_v,original_a,scaling_ok]=scaled_polynomials(c,h);
        if ~scaling_ok
            part.reason='polynomial_scaling_or_derivative_representation_unresolved';return;
        end
        % Batch only the original-coefficient endpoints. Normalized values
        % below are evaluated independently and retain the integrity check.
        batch_endpoints=isa(c,'double')&&isa(h,'double')&&~issparse(c)&&~issparse(h);
        if batch_endpoints
            original_q_ends=horner_endpoints(c,h);
            original_v_ends=horner_endpoints(original_v,h);
            original_a_ends=horner_endpoints(original_a,h);
        end
        qmin=Inf;qmax=-Inf;vmax=-Inf;amax=-Inf;
        previous_end=[];
        for j=1:pieces
            % Each physical derivative is scaled independently from the
            % ORIGINAL coefficients. Never form h^3 before multiplying by a,
            % or recover a derivative by dividing an underflowed q term.
            h0=h(j);p=P(j,:);vp=VP(j,:);ap=AP(j,:);
            uq=[0;1;stationary_roots(vp)];
            uv=[0;1;stationary_roots(ap)];
            qv=polyval(p,uq);vv=polyval(vp,uv);av=polyval(ap,[0;1]);
            if any(~isfinite([qv;vv;av])),part.reason='nonfinite_polynomial_evaluation';return;end
            % A second endpoint calculation uses original local-time
            % coefficients and Horner arithmetic, independently of scaling.
            if batch_endpoints
                original_ends=[original_q_ends(j,:)',original_v_ends(j,:)',original_a_ends(j,:)'];
            else
                original_ends=[polyval(c(j,:),[0;h0]),...
                    polyval(original_v(j,:),[0;h0]),polyval(original_a(j,:),[0;h0])];
            end
            % uq/uv start with 0 and 1, so these are the same Horner values
            % already computed above; stationary-root evaluation is unchanged.
            normalized_ends=[qv(1:2),vv(1:2),av];
            if any(~isfinite(original_ends(:)))
                part.reason='original_coefficient_endpoint_unresolved';return;
            end
            observed=max(abs(original_ends),[],1);
            discrepancy=max(abs(original_ends-normalized_ends),[],1);
            integrity_tol=[tol.position(axis),tol.velocity(axis),tol.acceleration(axis)]+1e-12*max(1,observed);
            endpoint_observed(axis,:)=max(endpoint_observed(axis,:),observed);
            endpoint_discrepancy(axis,:)=max(endpoint_discrepancy(axis,:),discrepancy);
            endpoint_tolerance(axis,:)=max(endpoint_tolerance(axis,:),integrity_tol);
            calculated=[max(abs(qv)),max(abs(vv)),max(abs(av))];
            if any(discrepancy>integrity_tol)||any(calculated+integrity_tol<observed)
                part.reason='independent_endpoint_check_unresolved';return;
            end
            % Keep the independently observed endpoints in the bound envelope
            % even when their difference is within the integrity tolerance.
            % Their matching local times preserve meaningful extremum times.
            qv=[qv;original_ends(:,1)];uq=[uq;0;1];
            vv=[vv;original_ends(:,2)];uv=[uv;0;1];
            av=[av;original_ends(:,3)];ua=[0;1;0;1];
            [lo,ix]=min(qv);if lo<qmin,qmin=lo;qtmin=breaks(j)+h0*uq(ix);end
            [hi,ix]=max(qv);if hi>qmax,qmax=hi;qtmax=breaks(j)+h0*uq(ix);end
            [hi,ix]=max(abs(vv));if hi>vmax,vmax=hi;vtmax=breaks(j)+h0*uv(ix);end
            [hi,ix]=max(abs(av));if hi>amax,amax=hi;atmax=breaks(j)+h0*ua(ix);end
            start_values=normalized_ends(1,:);
            end_values=normalized_ends(2,:);
            if ~isempty(previous_end)
                joinmax=max(joinmax,abs(start_values-previous_end));
            end
            previous_end=end_values;
        end
        extrema.q_min(axis)=qmin;extrema.q_max(axis)=qmax;
        extrema.velocity_abs_max(axis)=vmax;extrema.acceleration_abs_max(axis)=amax;
        when.q_min(axis)=qtmin;when.q_max(axis)=qtmax;
        when.velocity_abs_max(axis)=vtmax;when.acceleration_abs_max(axis)=atmax;
        q(:,axis)=reshape(ppval(seg.pp{axis},t),[],1);
    end
    if any(~isfinite(q(:))),part.reason='nonfinite_sample_evaluation';return;end
    joint_residual=zeros(1,2);
    if isfield(seg,'q_samples')
        if ~(isnumeric(seg.q_samples)&&isreal(seg.q_samples)&&isequal(size(seg.q_samples),[n 2])&&all(isfinite(seg.q_samples(:))))
            part.reason='invalid_stored_joint_samples';return;
        end
        joint_delta=abs(q-seg.q_samples);joint_residual=max(joint_delta,[],1);
    end
    q1=q(:,1);q2=q(:,2);
    if isempty(model) || model.is_default
        % Preserve the original verifier's rounding and output for 20/20.
        c20=cosd(20);s20=sind(20);
        u=c20*s20*(1+cosd(q2));v=s20*sind(q2);
        r=[c20*c20-s20*s20*cosd(q2),cosd(q1).*u-sind(q1).*v,sind(q1).*u+cosd(q1).*v];
    else
        % Independent FK: apply Rx(q1) Rz(a1) Rx(q2) Rz(a2) ex in
        % sequence, without calling the production FK or using cached
        % inverse-geometry constants. Angles and q are all in degrees.
        a=model.alpha_deg;c1=cosd(a(1));s1=sind(a(1));
        c2=cosd(a(2));s2=sind(a(2));
        after_q2=[repmat(c2,n,1),s2*cosd(q2),s2*sind(q2)];
        after_a1=[c1*after_q2(:,1)-s1*after_q2(:,2),...
            s1*after_q2(:,1)+c1*after_q2(:,2),after_q2(:,3)];
        r=[after_a1(:,1),cosd(q1).*after_a1(:,2)-sind(q1).*after_a1(:,3),...
            sind(q1).*after_a1(:,2)+cosd(q1).*after_a1(:,3)];
    end
    g=seg.gamma;target=[cosd(g(:,1)).*cosd(g(:,2)),cosd(g(:,1)).*sind(g(:,2)),-sind(g(:,1))];
    residual=max(sqrt(sum((r-target).^2,2)));
    margins=struct('joint_lower',extrema.q_min-bounds(:,1)',...
        'joint_upper',bounds(:,2)'-extrema.q_max,...
        'velocity',vlimit-extrema.velocity_abs_max,...
        'acceleration',alimit-extrema.acceleration_abs_max);
    excess=struct('joint_lower',max(-margins.joint_lower,0),...
        'joint_upper',max(-margins.joint_upper,0),...
        'velocity',max(-margins.velocity,0),'acceleration',max(-margins.acceleration,0));
    joint_ok=all(excess.joint_lower<=tol.position)&&all(excess.joint_upper<=tol.position);
    velocity_ok=all(excess.velocity<=tol.velocity);acceleration_ok=all(excess.acceleration<=tol.acceleration);
    % Relative guard is needed when tiny nonuniform intervals magnify roundoff.
    join_tolerance=[max(tol.position),max(tol.velocity),max(tol.acceleration)]+...
        1e-10*[max(abs([extrema.q_min extrema.q_max])),max(extrema.velocity_abs_max),max(extrema.acceleration_abs_max)];
    c2_ok=all(joinmax<=join_tolerance);
    part.extrema_resolved=true;part.extrema=extrema;part.extremum_time=when;
    part.margins=margins;part.raw_excess=excess;part.sample_direction_residual_max=residual;
    part.c2_join_max=joinmax;part.c2_join_tolerance=join_tolerance;
    part.endpoint_observed_abs_max=endpoint_observed;
    part.endpoint_discrepancy_max=endpoint_discrepancy;
    part.endpoint_consistency_tolerance=endpoint_tolerance;
    part.independent_endpoint_check_pass=true;
    part.sample_joint_residual_max=joint_residual;part.sample_joint_pass=all(joint_residual<=tol.sample_joint);
    part.position_bounds_pass=joint_ok;part.velocity_bounds_pass=velocity_ok;part.acceleration_bounds_pass=acceleration_ok;
    part.sample_fk_pass=isfinite(residual)&&residual<=tol.sample_direction_chord;
    part.c2_pass=c2_ok;
    part.passed=joint_ok&&velocity_ok&&acceleration_ok&&part.sample_fk_pass&&part.sample_joint_pass&&c2_ok;
    if part.passed,part.status='pass';else,part.status='violation';part.reason='bounds_fk_or_c2_check_failed';end
end

function model=verification_model(trajectory,options)
% Every available copy must describe the same model. A trajectory's model
% takes precedence over the default, never over a conflicting explicit one.
    model=[];
    if isfield(trajectory,'model')
        model=coneik_resolve_model(trajectory.model);
    end
    if isfield(trajectory,'options')&&isstruct(trajectory.options)&&...
            isscalar(trajectory.options)&&isfield(trajectory.options,'model')
        saved=coneik_resolve_model(trajectory.options.model);
        if ~isempty(model)&&~strcmp(saved.key,model.key)
            error('coneik_verify_trajectory:ModelMismatch','Stored trajectory and options models disagree.');
        end
        model=saved;
    end
    if isfield(options,'model')
        supplied=coneik_resolve_model(options.model);
        if ~isempty(model)&&~strcmp(supplied.key,model.key)
            error('coneik_verify_trajectory:ModelMismatch','Explicit model differs from the stored trajectory model.');
        end
        model=supplied;
    end
    % A legacy trajectory has no metadata and retains the original 20/20
    % contract. A key without a full model cannot be reconstructed safely.
    has_segment_key=isstruct(trajectory.segments)&&isfield(trajectory.segments,'model_key');
    if isempty(model)&&(isfield(trajectory,'model_key')||has_segment_key)
        error('coneik_verify_trajectory:ModelMismatch','A stored model key requires full model metadata or an explicit model.');
    end
    if ~isempty(model)
        if isfield(trajectory,'model_key')&&~isequal(trajectory.model_key,model.key)
            error('coneik_verify_trajectory:ModelMismatch','Stored trajectory model key is inconsistent.');
        end
        if has_segment_key
            for k=1:numel(trajectory.segments)
                if ~isequal(trajectory.segments(k).model_key,model.key)
                    error('coneik_verify_trajectory:ModelMismatch','A segment model key differs from the verification model.');
                end
            end
        end
    end
end

function ends=horner_endpoints(c,h)
    % Same multiplication/addition order as array polyval, one row per piece.
    % Do not use sum(c,2), powers of h, or scaled coefficients here.
    left=c(:,1);right=left;
    for k=2:size(c,2)
        left=0.*left+c(:,k);right=h.*right+c(:,k);
    end
    ends=[left,right];
end

function [P,VP,AP,original_v,original_a,ok]=scaled_polynomials(c,h)
    % Reject curves the unchanged public evaluator cannot represent safely.
    % This is conservative even when a cancellation could make values finite.
    original_v=[3*c(:,1),2*c(:,2),c(:,3)];
    original_a=[6*c(:,1),2*c(:,2)];
    P=zeros(size(c));VP=zeros(size(c,1),3);AP=zeros(size(c,1),2);
    ok=all(isfinite(original_v(:)))&&all(isfinite(original_a(:)));
    if ~ok,return;end
    for j=1:4
        [P(:,j),good]=scaled_product(c(:,j),h,4-j,1);ok=ok&&good;
    end
    for j=1:3
        [VP(:,j),good]=scaled_product(c(:,j),h,3-j,4-j);ok=ok&&good;
    end
    for j=1:2
        [AP(:,j),good]=scaled_product(c(:,j),h,2-j,(4-j)*(3-j));ok=ok&&good;
    end
end

function [y,ok]=scaled_product(x,h,power,factor)
    % Mantissas stay O(1); exponents are combined before the final rounding.
    % Nonzero terms that cannot be represented normally are unresolved,
    % never silently replaced by zero and allowed to pass a bound check.
    y=zeros(size(x));active=x~=0;
    if any(active)
        [fx,ex]=log2(abs(x(active)));[fh,eh]=log2(h(active));[ff,ef]=log2(factor);
        y(active)=sign(x(active)).*pow2(fx.*(fh.^power)*ff,ex+power*eh+ef);
    end
    ok=all(isfinite(y))&&~any(active&(y==0|abs(y)<realmin));
end

function r=stationary_roots(p)
    first=find(p~=0,1);
    if isempty(first)||first==numel(p),r=zeros(0,1);return;end
    p=p(first:end);original=p;p=p/max(abs(p));
    if any(original~=0&p==0)
        error('coneik_verify_trajectory:NumericalUnresolved','Stationary-root coefficient scaling lost a nonzero term.');
    end
    if numel(p)==2,r=-p(2)/p(1);
    else
        d=p(2)^2-4*p(1)*p(3);
        if d<0,r=zeros(0,1);return;end
        sd=sqrt(d);sgn=1;if p(2)<0,sgn=-1;end
        z=-0.5*(p(2)+sgn*sd);
        if z==0,r=-p(2)/(2*p(1));else,r=[z/p(1);p(3)/z];end
    end
    r=r(isfinite(r)&r>0&r<1);
end

function part=empty_part()
    part=struct('status','unresolved','passed',false,'reason','',...
        'extrema_resolved',false,'extrema',struct(),'extremum_time',struct(),...
        'margins',struct(),'raw_excess',struct(),'sample_direction_residual_max',NaN,...
        'c2_join_max',[NaN NaN NaN],'c2_join_tolerance',[NaN NaN NaN],...
        'endpoint_observed_abs_max',nan(2,3),'endpoint_discrepancy_max',nan(2,3),...
        'endpoint_consistency_tolerance',nan(2,3),'independent_endpoint_check_pass',false,...
        'sample_joint_residual_max',[NaN NaN],'sample_joint_pass',false,...
        'position_bounds_pass',false,'velocity_bounds_pass',false,...
        'acceleration_bounds_pass',false,'sample_fk_pass',false,'c2_pass',false);
end
function v=get_option(o,name,default)
    if isfield(o,name),v=o.(name);else,v=default;end
end
function validate_limits(b,v,a)
    if ~(isnumeric(b)&&isreal(b)&&isequal(size(b),[2 2])&&~any(isnan(b(:)))&&all(b(:,1)<=b(:,2))&&~any(b(:,1)==Inf)&&~any(b(:,2)==-Inf)&&...
            isnumeric(v)&&isreal(v)&&isvector(v)&&numel(v)==2&&all(v(:)>=0)&&...
            isnumeric(a)&&isreal(a)&&isvector(a)&&numel(a)==2&&all(a(:)>=0))
        error('coneik_verify_trajectory:InvalidLimits','Invalid joint, velocity or acceleration limits.');
    end
end

function cfg=tracking_options(o)
    cfg=struct('enabled',isfield(o,'intersample_direction_tolerance_deg'),...
        'degrees',NaN,'chord',NaN,'budget',1024);
    if ~cfg.enabled
        if isfield(o,'intersample_direction_max_refinements')
            error('coneik_verify_trajectory:InvalidOptions','A direction refinement budget requires an explicit tolerance.');
        end
        return;
    end
    v=o.intersample_direction_tolerance_deg;
    if ~(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v)&&v>0&&v<180)
        error('coneik_verify_trajectory:InvalidOptions','Direction tolerance must satisfy 0<degrees<180.');
    end
    cfg.degrees=double(v);cfg.chord=2*sind(cfg.degrees/2);
    if ~(isfinite(cfg.chord)&&cfg.chord>0)
        error('coneik_verify_trajectory:InvalidOptions','Direction tolerance is not representable.');
    end
    if isfield(o,'intersample_direction_max_refinements')
        v=o.intersample_direction_max_refinements;
        if ~(isnumeric(v)&&isreal(v)&&isscalar(v)&&isfinite(v)&&v>=0&&fix(v)==v&&v<=flintmax/2)
            error('coneik_verify_trajectory:InvalidOptions','Direction refinement budget must be a nonnegative integer.');
        end
        cfg.budget=double(v);
    end
end

function d=empty_tracking(cfg)
    d=struct('status','unresolved','reason','legacy_checks_not_passed',...
        'reference','shortest_arc_slerp_unit_direction_linear_time',...
        'tolerance_deg',cfg.degrees,'tolerance_chord',cfg.chord,...
        'max_refinements',cfg.budget,'refinements_used',0,...
        'reference_interval_count',0,'atomic_interval_count',0,...
        'evaluated_point_count',0,'observed_max_chord',0,'observed_max_angle_deg',0,...
        'model_upper_bound_chord',Inf,'model_upper_bound_angle_deg',180,...
        'worst_interval',[NaN NaN],'budget_exhausted',false,...
        'ambiguous_reference',false,'certified',false,...
        'bound_method','endpoint error + [rho*(A1+A2)+rho^2*(V1+V2)^2+omega_ref^2]*h^2/8 + float64 guards',...
        'partition_breaks',[],'max_numerical_guard_chord',0);
end

function d=verify_direction(seg,extrema,cfg,model)
% The existing verifier has already checked scalar cubics, original/scaled
% endpoints, derivative representation, stationary roots, sample FK and C2.
% Its segment-wide extrema also bound every atomic cell (conservatively).
    d=empty_tracking(cfg);d.reason='numerical_envelope_unresolved';rho=pi/180;
    t=double(seg.t(:));g=double(seg.gamma);
    target=[cosd(g(:,1)).*cosd(g(:,2)),cosd(g(:,1)).*sind(g(:,2)),-sind(g(:,1))];
    target=target./sqrt(sum(target.^2,2));
    d.reference_interval_count=numel(t)-1;
    breaks=cell(1,2);coeff=cell(1,2);cuts=t;
    for axis=1:2
        [breaks{axis},c,p,order]=unmkpp(seg.pp{axis});
        coeff{axis}=[zeros(p,4-order),double(c)];cuts=[cuts;breaks{axis}(:)]; %#ok<AGROW>
    end
    cuts=unique(cuts);d.partition_breaks=cuts';d.atomic_interval_count=numel(cuts)-1;
    V=extrema.velocity_abs_max;A=extrema.acceleration_abs_max;
    % Inflate only nonzero derivative bounds. Exact constants stay constant.
    V(V>0)=V(V>0)+64*eps(V(V>0));A(A>0)=A(A>0)+64*eps(A(A>0));
    [vterm,okv]=scaled_product(sum(V),rho,1,1);
    [vterm,okvv]=scaled_product(vterm,1,0,vterm);
    if sum(V)==0,vterm=0;okvv=true;end
    [aterm,oka]=scaled_product(sum(A),rho,1,1);
    if ~okv||~okvv||~oka||any(~isfinite([V A]))||~isfinite(vterm+aterm),return;end
    angles=[20 20];if ~isempty(model),angles=model.alpha_deg;end
    upper=0;worst=[NaN NaN];
    for j=1:numel(cuts)-1
        left=cuts(j);right=cuts(j+1);
        % Use the piece starting at/before LEFT, then evaluate that same
        % piece on both ends, including its LEFT limit at RIGHT.
        ix=[find(breaks{1}<=left,1,'last'),find(breaks{2}<=left,1,'last')];
        ri=find(t<=left,1,'last');h_ref=t(ri+1)-t(ri);
        a=target(ri,:);b=target(ri+1,:);cr=cross(a,b);sn=norm(cr);cs=max(-1,min(1,dot(a,b)));
        theta=atan2(sn,cs);approx=0;tangent=zeros(1,3);
        if sn<=sqrt(eps)
            if cs<0
                d.reason='ambiguous_antipodal_reference';d.ambiguous_reference=true;d.worst_interval=[left right];return;
            end
            approx=norm(b-a);ref_guard=64*eps;
        else
            tangent=cross(cr,a)/sn;tangent=tangent/norm(tangent);
            ref_guard=64*eps/sn;
        end
        omega=theta/h_ref;
        [wterm,okw]=scaled_product(omega,1,0,omega);
        if theta==0,wterm=0;okw=true;end
        if ~okw||~isfinite(omega)||(theta~=0&&(omega==0||abs(omega)<realmin)),return;end
        M=aterm+vterm+wterm;
        if ~isfinite(M),return;end
        if M>0,M=M+64*eps(M);end
        % Explicit stack gives finite memory/work without recursion limits.
        stack=[left right];
        while ~isempty(stack)
            cell_=stack(end,:);stack(end,:)=[];lo=cell_(1);hi=cell_(2);h=hi-lo;
            e=zeros(1,2);guard=e;
            for endpoint=1:2
                time=cell_(endpoint);q=zeros(1,2);qguard=q;
                for axis=1:2
                    local=time-breaks{axis}(ix(axis));c=coeff{axis}(ix(axis),:);
                    q(axis)=polyval(c,local);
                    magnitude=polyval(abs(c),abs(local));
                    qguard(axis)=64*eps(max(1,magnitude))+...
                        4*eps(max(abs([time breaks{axis}(ix(axis))])))*V(axis);
                end
                phase=(time-t(ri))/h_ref*theta;
                if sn<=sqrt(eps),ref=a;else,ref=cos(phase)*a+sin(phase)*tangent;end
                r=tracking_fk(q,angles);e(endpoint)=norm(r-ref);
                guard(endpoint)=rho*sum(qguard)+approx+ref_guard+64*eps;
                if any(~isfinite([q qguard phase ref r e(endpoint) guard(endpoint)]))
                    d.reason='nonfinite_direction_evaluation';d.worst_interval=cell_;return;
                end
                d.evaluated_point_count=d.evaluated_point_count+1;
            end
            d.observed_max_chord=max(d.observed_max_chord,max(e));
            d.observed_max_angle_deg=chord_angle(d.observed_max_chord);
            d.max_numerical_guard_chord=max(d.max_numerical_guard_chord,max(guard));
            if any(e>cfg.chord+guard)
                d.status='violation';d.reason='observed_direction_tolerance_exceeded';d.worst_interval=cell_;return;
            end
            [remainder,okr]=scaled_product(M,h,2,1/8);
            raw=max(e+guard)+remainder;bound=raw+64*eps(max(1,raw));
            if ~okr||~isfinite(bound)
                d.reason='numerical_envelope_unresolved';d.worst_interval=cell_;return;
            end
            if bound<=cfg.chord
                if bound>upper,upper=bound;worst=cell_;end
                continue;
            end
            if d.refinements_used>=cfg.budget
                d.reason='refinement_budget_exhausted';d.budget_exhausted=true;d.worst_interval=cell_;return;
            end
            mid=lo+h/2;
            if ~(isfinite(mid)&&mid>lo&&mid<hi)
                d.reason='midpoint_not_representable';d.worst_interval=cell_;return;
            end
            d.refinements_used=d.refinements_used+1;
            stack=[stack;mid hi;lo mid]; %#ok<AGROW>
        end
    end
    d.status='pass';d.reason='';d.model_upper_bound_chord=upper;
    d.model_upper_bound_angle_deg=chord_angle(upper);d.worst_interval=worst;
end

function r=tracking_fk(q,angles)
% Independent rotation-chain FK, without production inverse/FK helpers.
    ca=cosd(angles(1));sa=sind(angles(1));cb=cosd(angles(2));sb=sind(angles(2));
    y=sb*cosd(q(2));z=sb*sind(q(2));x=ca*cb-sa*y;y=sa*cb+ca*y;
    r=[x,cosd(q(1))*y-sind(q(1))*z,sind(q(1))*y+cosd(q(1))*z];
end

function a=chord_angle(c)
    a=2*asind(min(1,max(0,c)/2));
end

function d=aggregate_tracking(parts,cfg)
    d=empty_tracking(cfg);items=[parts.intersample_direction];
    names={'refinements_used','reference_interval_count','atomic_interval_count','evaluated_point_count'};
    for k=1:numel(names),d.(names{k})=sum([items.(names{k})]);end
    d.observed_max_chord=max([items.observed_max_chord]);d.observed_max_angle_deg=chord_angle(d.observed_max_chord);
    d.max_numerical_guard_chord=max([items.max_numerical_guard_chord]);
    d.budget_exhausted=any([items.budget_exhausted]);d.ambiguous_reference=any([items.ambiguous_reference]);
    statuses={items.status};index=find(strcmp(statuses,'unresolved'),1);
    if isempty(index),index=find(strcmp(statuses,'violation'),1);end
    if isempty(index)
        d.status='pass';d.reason='';
        [d.model_upper_bound_chord,index]=max([items.model_upper_bound_chord]);
        d.model_upper_bound_angle_deg=chord_angle(d.model_upper_bound_chord);
    else
        d.status=items(index).status;d.reason=items(index).reason;
    end
    d.worst_interval=items(index).worst_interval;
    % Partition breaks are per-segment only, to avoid suggesting gap coverage.
end
