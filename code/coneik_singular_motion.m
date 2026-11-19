function [v,a,report] = coneik_singular_motion(J,C,gd,ga,factor,v0,options,vtol,atol)
%CONEIK_SINGULAR_MOTION Internal rank-one joint velocity/acceleration search.
% v=v0+n*z; curvature and retained acceleration are quadratic in z.
% Intersect scalar quadratic inequalities for task residuals and box limits.
% The acceleration null direction is then used to satisfy its box. A provided
% current_velocity is fixed. Every result is checked in the original J/C.
% This is a float64 selected-configuration search, not a global certificate.
    v=nan(2,1);a=v;
    refining=isfield(options,'state_selection_context');best_score=Inf;
    report=struct('attempted',true,'method','rank-one polynomial interval search',...
        'status','unresolved','velocity_fixed',isfield(options,'current_velocity'),...
        'parameter_intervals',zeros(0,2),'tested_parameters',0,...
        'selected_parameter',NaN,'acceleration_null_parameter',NaN,...
        'certified',false,'global_feasibility_checked',false);
    if refining
        report.witness_count=0;report.witness_scores=zeros(1,0);
        report.selected_score=NaN;report.first_witness_retained=false;
    end
    if factor.rank~=1,return;end
    keep=find(factor.active); drop=find(~factor.active);
    u=factor.U(:,keep);t=factor.V(:,keep);sigma=factor.s(keep);
    n=factor.V(:,drop);[~,j]=max(abs(n));if n(j)<0,n=-n;end
    nv=n;if report.velocity_fixed,nv=zeros(2,1);end
    vl=limits(options,'velocity_limits');al=limits(options,'acceleration_limits');
    % Columns represent z^2, z, and 1. Do not sample a curvature polynomial.
    B=C*[nv(1)^2,2*v0(1)*nv(1),v0(1)^2;...
        nv(1)*nv(2),v0(1)*nv(2)+v0(2)*nv(1),v0(1)*v0(2);...
        nv(2)^2,2*v0(2)*nv(2),v0(2)^2];
    rhs=-B;rhs(:,3)=rhs(:,3)+ga;
    d=(u'*rhs)/sigma; A=t*d;
    residual=J*A+B;residual(:,3)=residual(:,3)-ga;
    inequalities=zeros(0,3);
    rv=J*[zeros(2,1),nv,v0];rv(:,3)=rv(:,3)-gd;
    for i=1:2
        inequalities=[inequalities;rv(i,:)-[0 0 vtol/2];-rv(i,:)-[0 0 vtol/2];...
            residual(i,:)-[0 0 atol/2];-residual(i,:)-[0 0 atol/2]]; %#ok<AGROW>
        if isfinite(vl(i))
            p=[0,nv(i),v0(i)];
            inequalities=[inequalities;p-[0 0 vl(i)];-p-[0 0 vl(i)]]; %#ok<AGROW>
        end
    end
    % Image of a symmetric acceleration box under the retained row t'.
    active=t~=0;
    % Reserve most of the existing 1e-7 box tolerance for final verification.
    capacity=sum(abs(t(active)).*(al(active)+2.5e-8));
    if isfinite(capacity)
        inequalities=[inequalities;d-[0 0 capacity];-d-[0 0 capacity]];
    end
    if ~all(isfinite([B(:);A(:);inequalities(:)])),numerical();end
    intervals=[-Inf,Inf];if report.velocity_fixed,intervals=[0,0];end
    for k=1:size(inequalities,1)
        intervals=intersect_sets(intervals,nonpositive(inequalities(k,:)));
        if isempty(intervals),break;end
    end
    report.parameter_intervals=intervals;
    if isempty(intervals),report.status='no_feasible_parameter';return;end
    % Prefer exact projected compatibility roots over tolerance-edge points.
    compatibility=factor.U(:,drop)'*rhs;
    roots_z=real_roots(compatibility);
    seeds=[0,reshape(roots_z,1,[])];
    for k=1:size(intervals,1)
        lo=intervals(k,1);hi=intervals(k,2);
        seeds(end+1)=max(lo,min(hi,0)); %#ok<AGROW>
        if isfinite(lo)&&isfinite(hi)
            seeds(end+1)=lo/2+hi/2; %#ok<AGROW>
        end
    end
    if refining && ~report.velocity_fixed
        preference=options.state_selection_context.preferred_velocity;
        if ~isempty(preference)
            preferred_z=nv'*(preference-v0);
            if isfinite(preferred_z)
                for k=1:size(intervals,1)
                    seeds(end+1)=max(intervals(k,1),min(intervals(k,2),preferred_z)); %#ok<AGROW>
                end
            end
        end
    end
    seeds=unique(seeds(isfinite(seeds)),'stable');
    for z=seeds
        if ~any(z>=intervals(:,1)&z<=intervals(:,2)),continue;end
        report.tested_parameters=report.tested_parameters+1;
        vv=v0+nv*z;
        if refining&&report.velocity_fixed,vv=v0;end % preserve measured bits, including -0
        bb=C*[vv(1)^2;vv(1)*vv(2);vv(2)^2];
        aa=t*((u'*(ga-bb))/sigma);
        try
            [w,ok,winterval]=box_parameter(aa,n,al);
        catch ex
            if ~refining||~strcmp(ex.identifier,'coneik_inverse_state:NumericalUnresolved'),rethrow(ex);end
            continue;
        end
        if ~ok,continue;end
        ws=w;
        if refining
            preference=options.state_selection_context.preferred_acceleration;
            if ~isempty(preference)
                preferred_w=n'*(preference-aa);
                if isfinite(preferred_w),ws(end+1)=max(winterval(1),min(winterval(2),preferred_w));end
            end
            if all(isfinite(winterval)),ws(end+1)=winterval(1)/2+winterval(2)/2;end
            ws=unique(ws,'stable');
        end
        retained_a=aa;
        for w=ws
            aa=retained_a+n*w;
            if ~all(isfinite([vv;aa;bb])),continue;end
            if norm(J*vv-gd,Inf)<=vtol && norm(J*aa+bb-ga,Inf)<=atol && ...
                    all(abs(vv)<=vl+1e-8) && all(abs(aa)<=al+1e-7)
                if refining
                    score=coneik_state_selection_score(options.state_selection_context,[],vv,aa);
                    report.witness_count=report.witness_count+1;
                    report.witness_scores(end+1)=score;
                    if score>=best_score,continue;end
                    best_score=score;report.selected_score=score;
                    report.first_witness_retained=report.witness_count==1;
                end
                v=vv;a=aa;report.status='feasible';report.selected_parameter=z;
                report.acceleration_null_parameter=w;
                if ~refining,return;end
            end
        end
    end
    if refining&&strcmp(report.status,'feasible'),return;end
    % Nonempty projected sets without a verified witness are unresolved.
    numerical();
end

function value=limits(options,name)
    value=[Inf;Inf];if isfield(options,name),value=options.(name);end
end

function [w,ok,interval]=box_parameter(a,n,limit)
    interval=[NaN,NaN];
    for slack=[0,2.5e-8]
        lo=-Inf;hi=Inf;
        for k=1:2
            if isinf(limit(k)),continue;end
            if n(k)==0
                if abs(a(k))>limit(k)+slack,lo=Inf;hi=-Inf;break;end
            else
                bound=limit(k)+slack;
                pair=sort([(-bound-a(k))/n(k),(bound-a(k))/n(k)]);
                if any(~isfinite(pair)),numerical();end
                lo=max(lo,pair(1));hi=min(hi,pair(2));
            end
        end
        if lo<=hi,w=max(lo,min(hi,0));ok=true;interval=[lo,hi];return;end
    end
    w=NaN;ok=false;
end

function out=intersect_sets(a,b)
    out=zeros(0,2);
    for i=1:size(a,1)
        for j=1:size(b,1)
            lo=max(a(i,1),b(j,1));hi=min(a(i,2),b(j,2));
            if lo<=hi,out(end+1,:)=[lo,hi];end %#ok<AGROW>
        end
    end
end

function intervals=nonpositive(p)
    if p(1)==0
        if p(2)==0
            intervals=zeros(0,2);if p(3)<=0,intervals=[-Inf Inf];end
        else
            boundary=-p(3)/p(2);if ~isfinite(boundary),numerical();end
            if p(2)>0,intervals=[-Inf boundary];else,intervals=[boundary Inf];end
        end
        return;
    end
    r=real_roots(p);
    if isempty(r)
        intervals=zeros(0,2);if p(1)<0,intervals=[-Inf Inf];end
    elseif numel(r)==1
        if p(1)>0,intervals=[r r];else,intervals=[-Inf Inf];end
    elseif p(1)>0,intervals=[r(1),r(2)];
    else,intervals=[-Inf,r(1);r(2),Inf];end
end

function r=real_roots(p)
    if ~all(isfinite(p)),numerical();end
    if p(1)==0
        r=[];if p(2)~=0,r=-p(3)/p(2);end
    else
        scale=max(abs(p));q=p/scale;
        if any(p~=0&q==0),numerical();end
        discriminant=q(2)^2-4*q(1)*q(3);
        % A small negative discriminant is not silently snapped to zero.
        if discriminant<0
            if abs(discriminant)<=32*eps*(q(2)^2+abs(4*q(1)*q(3))),numerical();end
            r=[];return;
        end
        if discriminant==0,r=-q(2)/(2*q(1));
        else
            sg=1;if q(2)<0,sg=-1;end
            h=-.5*(q(2)+sg*sqrt(discriminant));
            r=sort([h/q(1),q(3)/h]);
        end
    end
    if any(~isfinite(r)),numerical();end
end

function numerical()
    error('coneik_inverse_state:NumericalUnresolved',...
        'Rank-one joint motion search could not verify a finite feasible parameter.');
end
