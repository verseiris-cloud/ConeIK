function q=coneik_try_select_position_xyz_q(target,reference_q,c)
% Validated dense ordinary-tolerance XYZ distance query; both MEX available.
% Success-only native probe: no backend cache/discovery or precision solve.
% A decline returns NaN; the caller invokes the unchanged resolver once.
    q=nan(2,1);model=c.model;
    identity=(target(2)==0 && target(3)==0) || ...
        (target(1)==0 && model.outer_deg==90 && model.sum_roundoff_deg==0);
    if identity,return;end
    base=coneik_direction_double_mex(target,model.alpha_deg,c.tolerance);
    if ~all(isfinite(base)),return;end
    r=target/max(abs(target));r=r/norm(r);
    q=coneik_position_select_mex(base,r,reference_q,c.bounds,c.w,model.alpha_deg,c.tolerance);
end
