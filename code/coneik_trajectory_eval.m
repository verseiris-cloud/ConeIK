function [result, info] = coneik_trajectory_eval(trajectory, t_query)
%CONEIK_TRAJECTORY_EVAL Evaluate the stored joint curve and its derivatives.
% [result,info] = coneik_trajectory_eval(trajectory,t_query) returns rows
% [t,q1,q2,v1,v2,a1,a2], in seconds and nominal degrees. Queries outside
% accepted segment domains (including gaps) return NaN in columns 2:7.
% No extrapolation or interpolation across an unsuccessful segment occurs.
% Nonfinite in-domain q/v/a raises coneik_trajectory_eval:NumericalUnresolved.
    if ~(isstruct(trajectory) && isscalar(trajectory) && isfield(trajectory,'segments'))
        error('coneik_trajectory_eval:InvalidTrajectory','A trajectory with segments is required.');
    end
    if ~(isnumeric(t_query) && isreal(t_query) && (isvector(t_query)||isempty(t_query)) && all(isfinite(t_query(:))))
        error('coneik_trajectory_eval:InvalidTime','Query times must be a finite real vector.');
    end
    t_query=double(t_query(:));
    result=[t_query,nan(numel(t_query),6)];
    segment_index=zeros(numel(t_query),1);
    for k=1:numel(trajectory.segments)
        seg=trajectory.segments(k);
        if ~(isfield(seg,'t')&&isfield(seg,'pp')&&isnumeric(seg.t)&&isreal(seg.t)&&isvector(seg.t)) ||...
                numel(seg.t)<2 || any(~isfinite(seg.t)) || any(~isfinite(diff(seg.t))) || any(diff(seg.t)<=0) || numel(seg.pp)~=2
            error('coneik_trajectory_eval:InvalidTrajectory','Each accepted segment needs two or more samples and two polynomials.');
        end
        for axis=1:2
            [breaks,coefs,~,order,dim]=unmkpp(seg.pp{axis});
            if dim~=1 || order<1 || order>4 || ~isreal(coefs) || ~isreal(breaks) || any(~isfinite(coefs(:))) || any(~isfinite(breaks)) ||...
                    any(~isfinite(diff(breaks))) || any(diff(breaks)<=0) || breaks(1)~=seg.t(1) || breaks(end)~=seg.t(end)
                error('coneik_trajectory_eval:InvalidTrajectory','Finite polynomial domains must match their segment endpoints.');
            end
        end
        mask=t_query>=seg.t(1) & t_query<=seg.t(end);
        if any(segment_index(mask)~=0)
            error('coneik_trajectory_eval:OverlappingSegments','Accepted segment domains must not overlap.');
        end
        for axis=1:2
            pp=seg.pp{axis};
            result(mask,1+axis)=reshape(ppval(pp,t_query(mask)),[],1);
            result(mask,3+axis)=reshape(ppval(derivative(pp,1),t_query(mask)),[],1);
            result(mask,5+axis)=reshape(ppval(derivative(pp,2),t_query(mask)),[],1);
        end
        segment_index(mask)=k;
    end
    domains=zeros(numel(trajectory.segments),2);
    for k=1:numel(trajectory.segments),domains(k,:)=[trajectory.segments(k).t(1),trajectory.segments(k).t(end)];end
    domains=sortrows(domains,1);
    if any(domains(2:end,1)<=domains(1:end-1,2))
        error('coneik_trajectory_eval:OverlappingSegments','Accepted segment domains must not overlap.');
    end
    unresolved_rows=(segment_index>0) & any(~isfinite(result(:,2:7)),2);
    if any(unresolved_rows)
        error('coneik_trajectory_eval:NumericalUnresolved', ...
            'In-domain position, velocity, or acceleration evaluation is nonfinite.');
    end
    info=struct('segment_index',segment_index,'in_domain',segment_index>0,...
        'outside_or_gap_rows',find(segment_index==0),...
        'units','seconds; nominal degrees, degrees/s, degrees/s^2',...
        'derivatives_from_stored_coefficients',true);
end

function out=derivative(pp,n)
    [breaks,coefs,~,order,dim]=unmkpp(pp);
    if dim~=1 || order>4 || order<1 || any(~isfinite(coefs(:)))
        error('coneik_trajectory_eval:InvalidTrajectory','Expected finite scalar polynomials of degree at most three.');
    end
    for j=1:n
        if size(coefs,2)==1, coefs=zeros(size(coefs,1),1);
        else, coefs=coefs(:,1:end-1).*(size(coefs,2)-1:-1:1); end
    end
    out=mkpp(breaks,coefs);
end
