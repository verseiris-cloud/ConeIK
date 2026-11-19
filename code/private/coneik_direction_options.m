function [policy,tolerance]=coneik_direction_options(options)
% Shared options for public and prepared XYZ entry points.
    if ~(isstruct(options)&&isscalar(options)) || any(~ismember(fieldnames(options),{'on_error','direction_tolerance_deg'}))
        error('coneik_inverse_direction:InvalidOptions','Unknown option or non-scalar options structure.');
    end
    policy='report';tolerance=1e-8;
    if isfield(options,'on_error')
        v=options.on_error;
        if ~(ischar(v)&&isrow(v) || isstring(v)&&isscalar(v)) || ~any(strcmp(v,{'report','error'}))
            error('coneik_inverse_direction:InvalidOptions','on_error must be report or error.');
        end
        policy=char(v);
    end
    if isfield(options,'direction_tolerance_deg')
        tolerance=options.direction_tolerance_deg;
        if ~(isnumeric(tolerance)&&isreal(tolerance)&&isscalar(tolerance)&&isfinite(tolerance)&&tolerance>0&&tolerance<180)
            error('coneik_inverse_direction:InvalidOptions','direction_tolerance_deg must be finite and in (0,180).');
        end
        tolerance=double(tolerance);
    end
end
