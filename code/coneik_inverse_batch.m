function [Q, info] = coneik_inverse_batch(gamma, model, options)
%CONEIK_INVERSE_BATCH Checked bulk position IK with optional C++ acceleration.
% GAMMA is N-by-2 [pitch,yaw] in degrees, Q is N-by-4 [q11,q21,q12,q22].
% MODEL is from coneik_model; omitted/[] means 20/20. OPTIONS.backend is
% 'auto' (default), 'matlab' or 'cpp'. Auto uses the MEX if installed.
% OPTIONS.on_error is 'report' (default, failed rows NaN) or 'error'.
% C++ handles ordinary rows and guarded equal/unequal-cone axial neighbours;
% General difficult rows try native 50/80 then 100/200-digit stages; identities and
% undecided rows request MATLAB. backend=matlab never calls either MEX.
% The ORIGINAL target doubles and model are passed to the full solver,
% including its increasing-precision route; no approximate fallback exists.
% INFO.fallback means C++ delegated to MATLAB (including exact identities).
% INFO.used_precision records confirmed high-precision use on successful
% rows. Native 80/200-digit rows have route=cpp_precision and are counted in
% native_precision_count; native_fast_count counts binary64 only. Solver-failed rows have precision_known=false and precision_dps=NaN
% because the scalar solver throws before returning precision diagnostics.
% Post-validation rejection retains known pre-validation solver precision.
% This API solves independent positions; it does not select a motion path.
% With only Q requested, validation='none' avoids building diagnostic cells.
% OPTIONS.validation='direction' independently evaluates both returned
% branches. A row failing direction_tolerance_deg (default 1e-8 degrees)
% is rejected, Q becomes NaN, and status is residual_rejected. Unavailable
% floating-point assessment also rejects the row (fail-closed). This option
% does not relax geometric reachability or replace high-precision fallback.
% Positive native backend handles are reused after first resolution. Call
% coneik_refresh_backends after replacing MEX files or switching installations.
% INFO.assessment stores pre-rejection residuals and local sensitivity;
% INFO.solver_success stores the original solver flags. These extra fields
% are present only when validation='direction'. Diagnostic sensitivity and
% boundary margins are estimates, not interval proofs or probabilities.

    if nargin<2,model=[];end
    if nargin<3,options=struct();end
    if ~(isnumeric(gamma)&&isreal(gamma)&&ismatrix(gamma)&&~issparse(gamma)&&size(gamma,2)==2)
        error('coneik_inverse_batch:InvalidInput','gamma must be a full real N-by-2 numeric matrix.');
    end
    [backend,policy,validation,tolerance]=coneik_batch_options(options);
    model=coneik_resolve_model(model);gamma=double(gamma);
    if nargout>1
        [Q,info]=coneik_inverse_batch_resolved(gamma,model,backend,policy,validation,tolerance);
    else
        Q=coneik_inverse_batch_resolved(gamma,model,backend,policy,validation,tolerance);
    end
end
