function options=coneik_motion_options(options)
%CONEIK_MOTION_OPTIONS Internal validation shared by the two state adapters.
    if isfield(options,'singular_motion')
        mode=options.singular_motion;
        if isstring(mode)&&isscalar(mode)&&~ismissing(mode),mode=char(mode);end
        if ~(ischar(mode)&&isrow(mode)&&any(strcmp(mode,{'joint','minimum_norm'})))
            error('coneik_inverse_state:InvalidOptions','singular_motion must be joint or minimum_norm.');
        end
        options.singular_motion=mode;
    end
    for name={'velocity_limits','acceleration_limits'}
        key=name{1};if ~isfield(options,key),continue;end
        value=options.(key);
        if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&...
                ~any(isnan(value(:)))&&all(value(:)>=0))
            error('coneik_inverse_state:InvalidOptions','%s must be a nonnegative real 2-vector; Inf is permitted.',key);
        end
        options.(key)=double(value(:));
    end
end
