function solve=coneik_prepare_inverse(model,options)
%CONEIK_PREPARE_INVERSE Prepare repeated angle-direction queries for one model.
% SOLVE=CONEIK_PREPARE_INVERSE(MODEL,OPTIONS) validates and captures immutable
% geometry and the same options accepted by CONEIK_INVERSE_BATCH. Call
% Q=SOLVE(GAMMA) or [Q,INFO]=SOLVE(GAMMA) with N-by-2 pitch/yaw degrees.
% Both branches, diagnostics, validation and complete precision fallbacks
% share the ordinary batch engine. Every query still validates GAMMA.
% Rebuild SOLVE to change geometry/options; later edits of the caller's
% structs cannot change an existing solver. No targets/results are cached.
% Native handles are NOT captured: CONEIK_REFRESH_BACKENDS also affects
% existing prepared solvers. backend=matlab never resolves a native backend.
% Optional C++ availability is checked at query time, as with the batch API.
% Example: solve=coneik_prepare_inverse(coneik_model(15,25)); Q=solve([7,20]);
    if nargin<1,model=[];end
    if nargin<2,options=struct();end
    [backend,policy,validation,tolerance]=coneik_batch_options(options);
    model=coneik_resolve_model(model);
    solve=@query;
    function [Q,info]=query(gamma)
        if ~(isnumeric(gamma)&&isreal(gamma)&&ismatrix(gamma)&&~issparse(gamma)&&size(gamma,2)==2)
            error('coneik_inverse_batch:InvalidInput','gamma must be a full real N-by-2 numeric matrix.');
        end
        gamma=double(gamma);
        if nargout>1
            [Q,info]=coneik_inverse_batch_resolved(gamma,model,backend,policy,validation,tolerance);
        else
            Q=coneik_inverse_batch_resolved(gamma,model,backend,policy,validation,tolerance);
        end
    end
end
