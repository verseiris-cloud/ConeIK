function [Q,info]=coneik_inverse_batch_resolved(gamma,model,backend,policy,validation,tolerance)
% Private shared engine. Inputs and options were validated by the caller.
% Resolve backend handles through the refreshable registry on every query.
    n=size(gamma,1);
    available=false;
    if ~strcmp(backend,'matlab'),available=~isempty(coneik_native_backend('position'));end
    if strcmp(backend,'cpp')&&~available
        error('coneik_inverse_batch:NativeUnavailable','C++ MEX is not installed. Run build_coneik_mex or select backend=matlab/auto.');
    end
    use_cpp=available&&~strcmp(backend,'matlab');
    if nargout<2 && strcmp(validation,'none')
        Q=values_only(gamma,model,use_cpp,policy);
        return;
    end
    used='matlab';if use_cpp,used='cpp';end
    info=struct('model_key',model.key,'backend_requested',backend,'backend_used',used,...
        'status',{repmat({'numerical_unresolved'},n,1)},'success',false(n,1),...
        'fallback',false(n,1),'used_precision',false(n,1),'precision_known',false(n,1),'precision_dps',nan(n,1),...
        'error_id',{repmat({''},n,1)},'error_message',{repmat({''},n,1)},...
        'position_info',{cell(n,1)},'native_fast_count',0,'native_precision_count',0,'matlab_row_count',0);
    if use_cpp
        [Q,codes,native_precision,precision_unreachable,native_dps]=coneik_native_positions(gamma,model.alpha_deg);
        fast=codes==0 & ~native_precision; unreachable=codes==-1; invalid=codes==-3;
        info.success(fast)=true;info.status(fast)={'regular'};info.precision_dps(fast)=16;
        position=struct('status','regular','route','cpp_analytic',...
            'used_fallback',false,'precision_dps',16,'branch_count',2,'model_key',model.key);
        info.position_info(fast)={position};
        info.success(native_precision)=true;info.status(native_precision)={'regular'};
        info.used_precision(native_precision)=true;info.precision_dps(native_precision)=native_dps(native_precision);
        position.route='cpp_precision';position.used_fallback=true;
        if model.is_default,position=rmfield(position,'model_key');end
        for level=reshape(unique(native_dps(native_precision)),1,[])
            position.precision_dps=level;info.position_info(native_precision & native_dps==level)={position};
        end
        info.native_precision_count=nnz(native_precision);
        info.status(unreachable)={'unreachable'};
        info.error_id(unreachable)={'inverse_motion_analysis:Unreachable'};
        info.error_message(unreachable)={'Target direction is outside this model reachable annulus.'};
        info.error_message(precision_unreachable)={coneik_precision_unreachable_message(model.alpha_deg)};
        info.status(invalid)={'invalid_input'};
        info.error_id(invalid)={'inverse_motion_analysis:InvalidInput'};
        info.error_message(invalid)={'gamma1 and gamma2 must be finite real numeric scalars in degrees.'};
        info.fallback=codes==-2; rows=find(info.fallback);
        info.native_fast_count=nnz(fast);info.matlab_row_count=numel(rows);
    else
        Q=nan(n,4);rows=(1:n)';info.matlab_row_count=n;
    end
    for k=reshape(rows,1,[])
        try
            if model.is_default
                % The model is already validated for this batch. Retain
                % the original default call without repeated struct checks.
                [Q(k,1),Q(k,2),Q(k,3),Q(k,4),position]=...
                    coneik_inverse_default_resolved(gamma(k,1),gamma(k,2),false);
            else
                [Q(k,1),Q(k,2),Q(k,3),Q(k,4),position]=...
                    coneik_inverse_position_resolved(gamma(k,1),gamma(k,2),model,false);
            end
            info.success(k)=true;info.status{k}=position.status;
            info.used_precision(k)=position.used_fallback;
            info.precision_dps(k)=position.precision_dps;info.position_info{k}=position;
        catch exception
            switch exception.identifier
                case 'inverse_motion_analysis:Unreachable',status='unreachable';
                case 'inverse_motion_analysis:InvalidInput',status='invalid_input';
                case 'inverse_motion_analysis:NumericalUnresolved',status='numerical_unresolved';
                case 'inverse_motion_analysis:PrecisionUnavailable',status='precision_unavailable';
                otherwise, rethrow(exception);
            end
            Q(k,:)=NaN;info.status{k}=status;
            info.error_id{k}=exception.identifier;info.error_message{k}=exception.message;
        end
    end
    info.precision_known=info.success;
    if strcmp(validation,'direction')
        info.solver_success=info.success;
        info.assessment=coneik_assess_solution_resolved(gamma,Q,model,tolerance);
        rejected=info.success & ~all(info.assessment.direction_pass,2);
        Q(rejected,:)=NaN;
        info.success(rejected)=false;
        info.status(rejected)={'residual_rejected'};
        info.error_id(rejected)={'coneik_inverse_batch:ResidualRejected'};
        info.error_message(rejected)={'Independent direction residual exceeds the requested tolerance or cannot be evaluated.'};
    end
    if strcmp(policy,'error')
        failed=find(~info.success,1);
        if ~isempty(failed)
            error(info.error_id{failed},'Batch row %d: %s',failed,info.error_message{failed});
        end
    end
end

function Q=values_only(gamma,model,use_cpp,policy)
% Keep row ordering, strict classification and all precision fallbacks;
% only omit the unused public diagnostics. Still propagate unexpected errors.
    n=size(gamma,1);first_failed=Inf;failure_id='';failure_message='';
    if use_cpp
        [Q,codes,~,precision_unreachable]=coneik_native_positions(gamma,model.alpha_deg);
        rows=find(codes==-2);
        if strcmp(policy,'error')
            failed=find(codes==-1 | codes==-3,1);
            if ~isempty(failed)
                first_failed=failed;
                if codes(failed)==-1
                    failure_id='inverse_motion_analysis:Unreachable';
                    failure_message='Target direction is outside this model reachable annulus.';
                    if precision_unreachable(failed),failure_message=coneik_precision_unreachable_message(model.alpha_deg);end
                else
                    failure_id='inverse_motion_analysis:InvalidInput';
                    failure_message='gamma1 and gamma2 must be finite real numeric scalars in degrees.';
                end
            end
        end
    else
        Q=nan(n,4);rows=(1:n)';
    end
    for k=reshape(rows,1,[])
        try
            if model.is_default
                [Q(k,1),Q(k,2),Q(k,3),Q(k,4)]=coneik_inverse_default_resolved(gamma(k,1),gamma(k,2),false);
            else
                [Q(k,1),Q(k,2),Q(k,3),Q(k,4)]=coneik_inverse_position_resolved(gamma(k,1),gamma(k,2),model,false);
            end
        catch exception
            if ~any(strcmp(exception.identifier,{'inverse_motion_analysis:Unreachable','inverse_motion_analysis:InvalidInput',...
                    'inverse_motion_analysis:NumericalUnresolved','inverse_motion_analysis:PrecisionUnavailable'}))
                rethrow(exception);
            end
            Q(k,:)=NaN;
            if strcmp(policy,'error') && k<first_failed
                first_failed=k;failure_id=exception.identifier;failure_message=exception.message;
            end
        end
    end
    if isfinite(first_failed),error(failure_id,'Batch row %d: %s',first_failed,failure_message);end
end

