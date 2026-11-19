function model = coneik_resolve_model(model)
%CONEIK_RESOLVE_MODEL Resolve and validate an immutable CONEIK_MODEL value.
%   Omission or [] selects 20/20. All fields are checked against the
%   canonical constructor, including alpha, key and derived constants.
%   Persistent canonical values avoid repeating trigonometry for the same
%   geometry. A modified struct is rejected rather than trusted as a cache.

    persistent default_model last_model
    if isempty(default_model), default_model = coneik_model(); end
    if nargin == 0 || isempty(model)
        model = default_model;
        return;
    end
    if ~(isstruct(model) && isscalar(model) && isfield(model,'alpha_deg') && ...
            isnumeric(model.alpha_deg) && isreal(model.alpha_deg) && ...
            isequal(size(model.alpha_deg),[1,2]) && all(isfinite(model.alpha_deg)))
        error('coneik_resolve_model:InvalidModel', 'model must be an unmodified value returned by coneik_model.');
    end
    if isequal(model, default_model)
        model = default_model;
        return;
    elseif ~isempty(last_model) && isequal(model, last_model)
        model = last_model;
        return;
    end
    try
        canonical = coneik_model(model.alpha_deg(1),model.alpha_deg(2));
    catch cause
        error('coneik_resolve_model:InvalidModel', 'Invalid model geometry: %s',cause.message);
    end
    if ~isequal(model,canonical)
        error('coneik_resolve_model:InvalidModel', ...
            'Model fields were modified or use an unsupported schema; construct a fresh coneik_model.');
    end
    last_model = canonical;
    model = canonical;
end
