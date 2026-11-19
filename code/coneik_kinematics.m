function [r, gamma, J, b, detail] = coneik_kinematics(theta, dtheta, model)
%CONEIK_KINEMATICS Analytic kinematics, optionally using CONEIK_MODEL geometry.
%   An omitted/empty MODEL retains the original 20/20 arithmetic. The second
%   argument remains joint velocity; use [0;0] when only geometry is supplied.
%   Valid positive cone angles may still underflow in binary64. A geometry
%   with a zero computed cone sine raises NumericalUnresolved rather than
%   silently reporting a degenerate direction or zero motion.
%   THETA is a finite real numeric 2-vector of NOMINAL angles, in degrees.
%   Optional DTHETA has the same shape contract, in degrees/second; omission
%   means zero velocity. Inputs are converted to double. An explicitly empty
%   or non-2-vector input is invalid.
%   R is the unit output direction; GAMMA=[pitch;yaw] is its principal angle
%   representation in degrees; J=d(GAMMA)/d(THETA) is degree/degree.
%   B=Jdot*DTHETA is the curvature acceleration, in degrees/second^2:
%       gamma_dot = J*dtheta, gamma_ddot = J*ddtheta + b.
%   DETAIL exposes the direction Jacobian with respect to RADIAN angles and
%   compact 2x3 curvature coefficients. For v in degrees/second, curvature is
%       detail.curvature_coefficients * [v(1)^2; v(1)*v(2); v(2)^2].
%   Coefficients include all degree/radian conversions. No full Hessian
%   tensor is constructed. At most three outputs avoid second-order work.
%
%   Model: r=Rx(theta1)*Rz(20)*Rx(theta2)*Rz(20)*[1;0;0]. theta2 is already
%   nominal: do not add/subtract theta1 when using this derivative interface.

    if ~(isnumeric(theta) && isreal(theta) && isvector(theta) && ...
            numel(theta) == 2 && all(isfinite(theta(:))))
        error('coneik_kinematics:InvalidTheta', ...
            'theta must be a finite real numeric 2-vector in degrees.');
    end
    if nargin < 2
        dtheta = [0; 0];
    elseif ~(isnumeric(dtheta) && isreal(dtheta) && isvector(dtheta) && ...
            numel(dtheta) == 2 && all(isfinite(dtheta(:))))
        error('coneik_kinematics:InvalidVelocity', ...
            'dtheta must be a finite real numeric 2-vector in degrees/second.');
    end

    if nargin < 3, model=[]; end
    if ~isempty(model), model=coneik_resolve_model(model); end
    switch nargout
        case {0,1}, r=coneik_kinematics_resolved(theta,dtheta,model);
        case 2, [r,gamma]=coneik_kinematics_resolved(theta,dtheta,model);
        case 3, [r,gamma,J]=coneik_kinematics_resolved(theta,dtheta,model);
        case 4, [r,gamma,J,b]=coneik_kinematics_resolved(theta,dtheta,model);
        otherwise, [r,gamma,J,b,detail]=coneik_kinematics_resolved(theta,dtheta,model);
    end
end
