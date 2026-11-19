function [gamma1, gamma2, dgamma1, dgamma2] = theta_to_gamma(theta1, theta2, dtheta1, dtheta2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%THETA_TO_GAMMA Nominal angles/velocities to output angles/velocities.
%   Inputs are finite real numeric scalars in degrees and degrees/second.
%   Gamma and its degree/degree Jacobian share one analytic kernel call.

    if nargin < 5, model=[]; end
    values = {theta1,theta2,dtheta1,dtheta2};
    for k=1:numel(values)
        v=values{k};
        if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v))
            error('theta_to_gamma:InvalidInput', ...
                'Angles and velocities must be finite real numeric scalars.');
        end
    end
    theta=[double(theta1);double(theta2)];
    velocity=[double(dtheta1);double(dtheta2)];
    [~,gamma,J]=coneik_kinematics(theta,velocity,model);
    dgamma=J*velocity;
    gamma1=gamma(1); gamma2=gamma(2);
    dgamma1=dgamma(1); dgamma2=dgamma(2);
end
