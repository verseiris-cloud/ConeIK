function [backend,policy,validation,tolerance]=coneik_batch_options(options)
% Shared public option validation; never resolve or retain a native handle.
    if ~(isstruct(options)&&isscalar(options)) || any(~ismember(fieldnames(options),{'backend','on_error','validation','direction_tolerance_deg'}))
        error('coneik_inverse_batch:InvalidOptions','Unknown option or non-scalar options structure.');
    end
    backend=text_option(options,'backend','auto',{'auto','matlab','cpp'});
    policy=text_option(options,'on_error','report',{'report','error'});
    validation='none';
    if isfield(options,'validation'),validation=text_option(options,'validation','none',{'none','direction'});end
    tolerance=1e-8;
    if isfield(options,'direction_tolerance_deg')
        tolerance=options.direction_tolerance_deg;
        if ~strcmp(validation,'direction') || ~(isnumeric(tolerance)&&isreal(tolerance)&&isscalar(tolerance)&&isfinite(tolerance)&&tolerance>0&&tolerance<180)
            error('coneik_inverse_batch:InvalidOptions','direction_tolerance_deg requires validation=direction and a finite scalar in (0,180) degrees.');
        end
        tolerance=double(tolerance);
    end
end

function value=text_option(options,name,default,allowed)
    value=default;if isfield(options,name),value=options.(name);end
    if isstring(value)&&isscalar(value)&&~ismissing(value),value=char(value);end
    if ~(ischar(value)&&isrow(value)&&any(strcmp(value,allowed)))
        error('coneik_inverse_batch:InvalidOptions','Invalid %s option.',name);
    end
end
