function [q11,q21,q12,q22,info] = coneik_inverse_position(gamma1,gamma2,model)
%CONEIK_INVERSE_POSITION Parameterized stable half-angle direction inverse.
%   [Q11,Q21,Q12,Q22,INFO]=CONEIK_INVERSE_POSITION(PITCH,YAW,MODEL).
%   All inputs/outputs are degrees; Q2 is NOMINAL. The negative Q2 branch
%   comes first before periodic representative adjustment. Q1 is wrapped to
%   [-180,180); Q1+Q2 is represented in [-180,180]. MODEL is from CONEIK_MODEL.
%   Omission/[] preserves the original 20/20 solver exactly.
%   INFO separates regular/center/inner_boundary/outer_boundary and records
%   branch_count, route, used_fallback, precision_dps and model_key.
%   Invalid, unreachable and numerically unresolved inputs raise distinct
%   errors. No neighborhood is snapped onto a boundary. Difficult binary64
%   targets AND geometry are recomputed at increasing symbolic precision.

    if nargin < 3, model = []; end
    model = coneik_resolve_model(model);
    if model.is_default
        if nargout>4
            [q11,q21,q12,q22,info] = coneik_inverse_default_resolved(gamma1,gamma2,true);
            info.model_key = model.key;
        else
            [q11,q21,q12,q22] = coneik_inverse_default_resolved(gamma1,gamma2,true);
        end
        return;
    end
    if nargout>4
        [q11,q21,q12,q22,info]=coneik_inverse_position_resolved(gamma1,gamma2,model);
    else
        [q11,q21,q12,q22]=coneik_inverse_position_resolved(gamma1,gamma2,model);
    end
end
