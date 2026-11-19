function [Q,info]=coneik_inverse_batch_gpu(gamma,model,options)
%CONEIK_INVERSE_BATCH_GPU GPU batch twin of coneik_inverse_batch.
% [Q,INFO]=CONEIK_INVERSE_BATCH_GPU(GAMMA,MODEL,OPTIONS) converts the N-by-2
% [pitch,yaw] degree matrix to XYZ directions and evaluates them with the
% GPU half-angle core. Rows the GPU core declines (exact identities,
% boundaries, guarded or difficult inputs) fall back to the unchanged
% CONEIK_INVERSE_BATCH row engine, so outputs agree with the CPU API for
% every input. OPTIONS are the same as CONEIK_INVERSE_BATCH and are passed
% through to the fallback call; the backend setting is irrelevant on GPU
% rows. Requires the Parallel Computing Toolbox and a usable GPU device.
% Q is N-by-4 [q11,q21,q12,q22] in nominal degrees; INFO is the reduced
% diagnostic structure of the GPU path. This API solves independent
% positions; it does not select a motion path.
% Example:
%   Q=coneik_inverse_batch_gpu([20 -5; 0 0],coneik_model(20,20));
    if nargin<2,model=[];end
    if nargin<3,options=struct();end
    coneik_assert_gpu();
    if ~(isnumeric(gamma)&&isreal(gamma)&&ismatrix(gamma)&&~issparse(gamma)&&size(gamma,2)==2)
        error('coneik_inverse_batch:InvalidInput','gamma must be a full real N-by-2 numeric matrix.');
    end
    % The GPU gate uses the direction tolerance; CONEIK_INVERSE_BATCH fully
    % validates the complete options structure on the fallback call.
    tolerance=1e-8;
    if isstruct(options)&&isscalar(options)&&isfield(options,'direction_tolerance_deg')
        t=options.direction_tolerance_deg;
        if ~(isnumeric(t)&&isreal(t)&&isscalar(t)&&isfinite(t)&&t>0&&t<180)
            error('coneik_inverse_direction:InvalidOptions','direction_tolerance_deg must be finite and in (0,180).');
        end
        tolerance=double(t);
    end
    model=coneik_resolve_model(model);
    gamma=double(gamma);
    n=size(gamma,1);Q=nan(n,4);fast=false(n,1);residual=nan(n,2);
    if n>0
        % Same direction construction as the scalar position engine.
        cp=cosd(gamma(:,1));
        R=[cp.*cosd(gamma(:,2)),cp.*sind(gamma(:,2)),-sind(gamma(:,1))];
        [fastQ,fast,fastResidual]=coneik_direction_gpu_batch(R,model,tolerance);
        Q(fast,:)=fastQ(fast,:);residual(fast,:)=fastResidual(fast,:);
        if any(~fast)
            rows=find(~fast);
            % Unchanged row engine: exact identities, boundaries, axis
            % neighbours, precision ladder and error policy are all CPU.
            Q(rows,:)=coneik_inverse_batch(gamma(rows,:),model,options);
        end
    end
    if nargout>1
        ok=all(isfinite(Q),2);
        status=repmat({'regular'},n,1);status(~ok)={'cpu_row_engine_failed'};
        route=repmat({'vector_gpu_half_angle'},n,1);route(~fast)={'cpu_row_engine'};
        info=struct('model_key',model.key,'input_representation','pitch_yaw',...
            'success',ok,'solver_success',ok,'status',{status},...
            'route',{route},'precision_dps',16*ones(n,1),'used_precision',false(n,1),...
            'branch_count',2*ones(n,1),'direction_error_deg',residual,...
            'error_id',{repmat({''},n,1)},'error_message',{repmat({''},n,1)});
    end
end
