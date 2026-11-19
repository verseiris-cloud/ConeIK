function [Q, meta] = coneik_position_candidates(gamma, options, cache)
%CONEIK_POSITION_CANDIDATES Analytic branches and bounded periodic representatives.
% Q is 2-by-K in nominal degrees. OPTIONS accepts previous_q, joint_bounds
% (2-by-2, one [lower upper] row per joint), center_q1 (finite seeds),
% max_candidates (default 256) and winding_radius (default 1).
% All lifts inside two finite bounds are enumerated. An unbounded or
% one-sided axis uses winding_radius turns around the nearest admissible
% representative. Center q1 is a continuous family: its seeds are incomplete.
% Optional CACHE must be an unmodified META.cache from this function for the
% identical GAMMA and MODEL; it avoids repeating analytic IK during a path search.
% model is an optional coneik_model structure (default 20/20 degrees).
% Position bounds use the reported 1e-9 degree engineering tolerance.

    gamma=vector2(gamma,'gamma');
    if nargin<2, options=struct(); end
    allowed={'previous_q','joint_bounds','center_q1','max_candidates','winding_radius','model'};
    if ~(isstruct(options)&&isscalar(options)) || any(~ismember(fieldnames(options),allowed))
        error('coneik_position_candidates:InvalidOptions','Unknown option or non-scalar options struct.');
    end
    model=[];
    if isfield(options,'model')
        model=coneik_resolve_model(options.model);
    end
    if nargin<3,cache=[];end
    [Q,meta]=coneik_position_candidates_resolved(gamma,options,cache,model);
end

function value=vector2(value,name)
    if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&all(isfinite(value(:))))
        error('coneik_position_candidates:InvalidInput','%s must be a finite real 2-vector.',name);
    end
    value=double(value(:));
end
