function c=coneik_position_selection_options(model,options)
    allowed={'joint_bounds','weights','objective','on_error','direction_tolerance_deg'};
    if ~(isstruct(options)&&isscalar(options))||any(~ismember(fieldnames(options),allowed))
        error('coneik_select_inverse_position:InvalidOptions','Unknown option or non-scalar options.');
    end
    bounds=[-Inf,Inf;-Inf,Inf];if isfield(options,'joint_bounds'),bounds=options.joint_bounds;end
    if ~(isnumeric(bounds)&&isreal(bounds)&&isequal(size(bounds),[2,2])&&~any(isnan(bounds(:)))&&...
            all(bounds(:,1)<=bounds(:,2))&&~any(bounds(:,1)==Inf)&&~any(bounds(:,2)==-Inf))
        error('coneik_select_inverse_position:InvalidOptions','joint_bounds must contain valid [lower upper] rows.');
    end
    weights=[1;1];if isfield(options,'weights'),weights=options.weights;end
    if ~(isnumeric(weights)&&isreal(weights)&&isvector(weights)&&numel(weights)==2&&all(isfinite(weights(:)))&&...
            all(weights(:)>=0)&&any(weights(:)>0))
        error('coneik_select_inverse_position:InvalidOptions','weights must be finite, nonnegative and not both zero.');
    end
    weights=double(weights(:));weight_scale=max(weights);w=weights/weight_scale;
    if any(weights>0 & w==0)
        error('coneik_select_inverse_position:InvalidOptions','Weight ratio exceeds reliable double resolution.');
    end
    objective=text_option(options,'objective','distance',{'distance','ikjl_periodic'});
    policy=text_option(options,'on_error','report',{'report','error'});
    tolerance=1e-8;if isfield(options,'direction_tolerance_deg'),tolerance=options.direction_tolerance_deg;end
    if ~(isnumeric(tolerance)&&isreal(tolerance)&&isscalar(tolerance)&&isfinite(tolerance)&&tolerance>0&&tolerance<180)
        error('coneik_select_inverse_position:InvalidOptions','direction_tolerance_deg must be finite and in (0,180).');
    end
    bounds=double(bounds);tolerance=double(tolerance);
    model=coneik_resolve_model(model);
    c=struct('model',model,'bounds',bounds,'weights',weights,'weight_scale',weight_scale,'w',w,...
        'objective',objective,'policy',policy,'tolerance',tolerance);
end
function value=text_option(options,key,default,choices)
    value=default;if isfield(options,key),value=options.(key);end
    if isstring(value)&&isscalar(value)&&~ismissing(value),value=char(value);end
    if ~(ischar(value)&&isrow(value)&&any(strcmp(value,choices)))
        error('coneik_select_inverse_position:InvalidOptions','Invalid %s.',key);
    end
end
