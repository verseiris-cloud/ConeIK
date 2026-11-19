function assessment=coneik_assess_solution(gamma,Q,model,direction_tolerance_deg)
%CONEIK_ASSESS_SOLUTION Independently assess direction IK candidates.
% GAMMA is N-by-2 [pitch,yaw] degrees; Q is N-by-2 or N-by-4 NOMINAL
% joint angles in degrees. Nonfinite rows yield unavailable diagnostics and
% direction_pass=false. MODEL is an unmodified CONEIK_MODEL (or []/omitted).
% Tolerance defaults to 1e-8 degrees and must lie strictly in (0,180).
% Does not call the inverse solver, select branches, clamp targets or prove
% reachability. Its forward rotation chain is independent of the inverse.
%
% Fields with one column per branch:
%   direction_error_deg / direction_pass: FK residual and tolerance test.
%   unit_norm_error: absolute deviation of computed FK norm from one.
%   sigma_min / sigma_max: singular values of d(direction)/d(q_rad).
%   condition_number: sigma_max/sigma_min, Inf at numerical rank loss.
%   inverse_gain: 1/sigma_min, local joint-vector error / angular direction
%     error (rad/rad or deg/deg); not a finite-perturbation bound.
%   numerically_singular: sigma_min <= 32*eps*sigma_max, diagnostic only.
%   evaluation_available: finite inputs and representable direction geometry.
% boundary_margin_deg is N-by-1: approximate min distance to inner/outer
% cones. Its sign is NOT certified near a boundary. It never overrides
% the inverse solver's strict classification. No confidence probability or
% interval certificate is supplied, and sensor/model uncertainty is absent.
% A small residual does not imply unique or noise-insensitive joint angles.

    if nargin<3,model=[];end
    if nargin<4,direction_tolerance_deg=1e-8;end
    if ~(isnumeric(gamma)&&isreal(gamma)&&ismatrix(gamma)&&~issparse(gamma)&&size(gamma,2)==2)
        error('coneik_assess_solution:InvalidInput','gamma must be a full real N-by-2 numeric matrix.');
    end
    if ~(isnumeric(Q)&&isreal(Q)&&ismatrix(Q)&&~issparse(Q)&&size(Q,1)==size(gamma,1)&&any(size(Q,2)==[2 4]))
        error('coneik_assess_solution:InvalidInput','Q must be a full real N-by-2 or N-by-4 matrix matching gamma rows.');
    end
    if ~(isnumeric(direction_tolerance_deg)&&isreal(direction_tolerance_deg)&&isscalar(direction_tolerance_deg)&&...
            isfinite(direction_tolerance_deg)&&direction_tolerance_deg>0&&direction_tolerance_deg<180)
        error('coneik_assess_solution:InvalidTolerance','Tolerance must be a finite scalar in (0,180) degrees.');
    end
    model=coneik_resolve_model(model);
    assessment=coneik_assess_solution_resolved(double(gamma),double(Q),model,direction_tolerance_deg);
end
