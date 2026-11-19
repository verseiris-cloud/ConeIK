function solve=coneik_prepare_direction_gpu(model,options)
%CONEIK_PREPARE_DIRECTION_GPU GPU batch twin of coneik_prepare_direction.
% SOLVE=CONEIK_PREPARE_DIRECTION_GPU(MODEL,OPTIONS) captures the validated
% model and the same options as CONEIK_INVERSE_DIRECTION, then evaluates
% Q=SOLVE(R) for a full real N-by-3 direction matrix on the GPU. Ordinary
% interior rows run the vectorized half-angle core as gpuArray element
% operations. Every row the GPU core declines (identities, boundaries,
% ill-conditioned or tight-tolerance inputs) falls back to the unchanged
% CPU row engine, so outputs agree with CONEIK_PREPARE_DIRECTION for every
% input. Requires the Parallel Computing Toolbox and a usable GPU device.
% Angles are nominal degrees; Q contains both branches [q11,q21,q12,q22].
% Requesting INFO returns reduced diagnostics: route vector_gpu_half_angle
% marks GPU rows, cpu_row_engine marks fallback rows. There is no cache.
% Example:
%   s=coneik_prepare_direction_gpu(coneik_model(20,20));
%   Q=s([1,.25,-.125;1,-.2,.1]);
    if nargin<1,model=[];end
    if nargin<2,options=struct();end
    coneik_assert_gpu();
    [policy,tolerance]=coneik_direction_options(options);
    model=coneik_resolve_model(model);
    solve=@query;
    function [Q,info]=query(R)
        if ~(isnumeric(R)&&isreal(R)&&ismatrix(R)&&~issparse(R)&&size(R,2)==3)
            error('coneik_inverse_direction:InvalidInput','R must be a full real N-by-3 numeric matrix.');
        end
        R=double(R);
        n=size(R,1);Q=nan(n,4);fast=false(n,1);residual=nan(n,2);
        if n>0
            [fastQ,fast,fastResidual]=coneik_direction_gpu_batch(R,model,tolerance);
            Q(fast,:)=fastQ(fast,:);residual(fast,:)=fastResidual(fast,:);
        end
        if any(~fast)
            rows=find(~fast);
            % Unchanged CPU row engine keeps identities, boundaries and the
            % precision ladder, and applies the same error policy.
            Q(rows,:)=coneik_direction_resolved(R(rows,:),model,policy,tolerance);
        end
        if nargout>1
            ok=all(isfinite(Q),2);
            status=repmat({'regular'},n,1);status(~ok)={'cpu_row_engine_failed'};
            route=repmat({'vector_gpu_half_angle'},n,1);route(~fast)={'cpu_row_engine'};
            info=struct('model_key',model.key,'input_representation','direction_xyz',...
                'success',ok,'solver_success',ok,'status',{status},...
                'route',{route},'precision_dps',16*ones(n,1),'used_precision',false(n,1),...
                'branch_count',2*ones(n,1),'direction_error_deg',residual,...
                'error_id',{repmat({''},n,1)},'error_message',{repmat({''},n,1)});
        end
    end
end
