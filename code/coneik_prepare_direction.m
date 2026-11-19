function solve=coneik_prepare_direction(model,options)
%CONEIK_PREPARE_DIRECTION Reuse one geometry for repeated original-XYZ queries.
% SOLVE=CONEIK_PREPARE_DIRECTION(MODEL,OPTIONS) captures the validated model
% and the same options as CONEIK_INVERSE_DIRECTION. Call Q=SOLVE(R) or
% [Q,INFO]=SOLVE(R) with a real full N-by-3 direction matrix. Angles are
% nominal degrees; Q contains both branches [q11,q21,q12,q22].
% Both branches always undergo direction validation, including Q-only calls.
% Q-only calls avoid constructing unused public diagnostics; requesting INFO
% retains the complete diagnostic path and its original precision provenance.
% Ordinary scalar targets may use the optional vector_cpp_double route;
% absence or explicit decline retains the original MATLAB/fallback chain.
% Targets are validated on every query; difficult rows retain the original
% XYZ high-precision fallback, first trying the optional native 50/80 then 100/200/400
% backend. Exact 45-degree boundaries may use its exact dyadic predicate;
% unhandled identities/unconfirmed rows retain Symbolic fallback. This closure
% captures no MEX handle: CONEIK_REFRESH_BACKENDS affects subsequent queries.
% There is no target/result cache.
% Model/options are immutable copies. Rebuild SOLVE to change either.
% Example: solve=coneik_prepare_direction(coneik_model(15,25)); Q=solve([1,.25,-.125]);
% Optional native small-batch acceleration applies only to 8..16 XYZ rows
% at direction_tolerance_deg >= 1e-10; route vector_cpp_double_batch.
% Both branches retain residual checks. Declined rows retry the original
% MATLAB batch solver before its original precision/error fallback.
% Other sizes keep their established paths; missing MEX remains supported.
% CONEIK_REFRESH_BACKENDS is for backend/path changes, not per-target use.
    if nargin<1,model=[];end
    if nargin<2,options=struct();end
    [policy,tolerance]=coneik_direction_options(options);
    model=coneik_resolve_model(model);
    solve=@query;
    function [Q,info]=query(R)
        if ~(isnumeric(R)&&isreal(R)&&ismatrix(R)&&~issparse(R)&&size(R,2)==3)
            error('coneik_inverse_direction:InvalidInput','R must be a full real N-by-3 numeric matrix.');
        end
        R=double(R);
        if nargout>1
            [Q,info]=coneik_direction_resolved(R,model,policy,tolerance);
        else
            Q=coneik_direction_resolved(R,model,policy,tolerance);
        end
    end
end
