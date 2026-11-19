function [theta1, theta2, dtheta1, dtheta2] = gamma_to_theta(gamma1, gamma2, dgamma1, dgamma2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%GAMMA_TO_THETA First-branch inverse of direction angles and velocity.
%   Scalar angles/velocities are in degrees and degrees/second. Position is
%   delegated to the unchanged stable analytic inverse. The legacy unique-
%   inverse rejection (abs(det(J))<1e-12) is retained; singular projection
%   is not substituted. Nonprincipal pitch derivatives are converted into
%   the principal chart used by the forward Jacobian.

    if nargin < 5, model=[]; end
    [theta1,theta2,~,~]=inverse_motion_analysis(gamma1,gamma2,model);
    if cosd(rem(double(gamma1),360))==0
        error('gamma_to_theta:ChartSingularity','Pitch/yaw velocity is undefined at an exact direction pole.');
    end
    if ~(isnumeric(dgamma1) && isreal(dgamma1) && isscalar(dgamma1) && ...
            isfinite(dgamma1) && isnumeric(dgamma2) && isreal(dgamma2) && ...
            isscalar(dgamma2) && isfinite(dgamma2))
        error('gamma_to_theta:InvalidVelocity', ...
            'Direction velocities must be finite real numeric scalars.');
    end
    J=jacobian_gamma_theta(theta1,theta2,model);
    state=coneik_jacobian_state(J);
    if state.rank<2 || abs(det(J))<1e-12
        error('gamma_to_theta:NonUniqueInverse', ...
            'The configuration is singular or too close to singular for a unique velocity inverse.');
    end
    target=[double(dgamma1);double(dgamma2)];
    if cosd(rem(double(gamma1),360))<0
        target(1)=-target(1);
    end
    dtheta=J\target;
    tolerance=1e-8+1e-10*norm(target,inf);
    residual=J*dtheta-target;
    if ~all(isfinite(dtheta)) || ~all(isfinite(residual)) || norm(residual,inf)>tolerance
        error('gamma_to_theta:NumericalBreakdown', ...
            'The full-rank velocity solve did not meet the task residual tolerance.');
    end
    dtheta1=dtheta(1); dtheta2=dtheta(2);
end
