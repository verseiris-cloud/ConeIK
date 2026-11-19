function [Q,info]=coneik_inverse_direction(R,model,options)
%CONEIK_INVERSE_DIRECTION Invert N-by-3 beam directions without an angle chart.
% R contains finite, nonzero [x,y,z] rows; positive scale is immaterial.
% Q contains NOMINAL degree angles [q11,q21,q12,q22], negative-q2 branch first
% before periodic representative adjustment. The coordinate convention is
% R ~ [cosd(pitch)*cosd(yaw), cosd(pitch)*sind(yaw), -sind(pitch)].
% Original binary64 components define the target. Difficult rows recompute
% directly from these components, never from rounded pitch/yaw values.
% MODEL is from coneik_model. OPTIONS.on_error is 'report' (default, NaN rows)
% or 'error'. OPTIONS.direction_tolerance_deg defaults to 1e-8. Both returned
% branches are always checked against the input direction. Validation measures
% floating-point model closure, not measurement/geometry uncertainty.
% Q-only calls omit unused public diagnostics; both-branch validation,
% reachability decisions and all precision fallbacks still run unchanged.
% Ordinary single rows try the optional C++ double backend (vector_cpp_double)
% before the MATLAB double path. Batches retain the MATLAB array path.
% Difficult rows first try the
% optional native 50/80 then 100/200/400-digit backend. Confirmed precision_dps
% is 80, 200 or 400; skipped/unconfirmed low tiers retain the higher ladder.
% Geometry constants are reused; original XYZ is retained. Exact 45-degree
% inner/outer boundaries with an exact raw-XYZ Pythagorean relation can use
% vector_cpp_exact_boundary (80-digit phase recovery, exact classification).
% Other exact boundaries retain Symbolic fallback and may require its toolbox.
% Optional native small-batch acceleration applies only to 8..16 XYZ rows
% at direction_tolerance_deg >= 1e-10; route vector_cpp_double_batch.
% Both branches retain residual checks. Declined rows retry the original
% MATLAB batch solver before its original precision/error fallback.
% Other sizes keep their established paths; missing MEX remains supported.
% CONEIK_REFRESH_BACKENDS is for backend/path changes, not per-target use.
    if nargin<2,model=[];end
    if nargin<3,options=struct();end
    if ~(isnumeric(R)&&isreal(R)&&ismatrix(R)&&~issparse(R)&&size(R,2)==3)
        error('coneik_inverse_direction:InvalidInput','R must be a full real N-by-3 numeric matrix.');
    end
    [policy,tolerance]=coneik_direction_options(options);
    model=coneik_resolve_model(model);R=double(R);
    if nargout>1
        [Q,info]=coneik_direction_resolved(R,model,policy,tolerance);
    else
        Q=coneik_direction_resolved(R,model,policy,tolerance);
    end
end
