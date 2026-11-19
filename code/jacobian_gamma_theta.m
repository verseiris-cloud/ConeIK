function J = jacobian_gamma_theta(theta1, theta2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%JACOBIAN_GAMMA_THETA Analytic d(pitch,yaw)/d(nominal theta1,theta2).
%   Inputs are finite real numeric scalars in degrees. J is degree/degree;
%   both input and output coordinate conversions cancel in this first
%   derivative. Curvature requires the separate coneik_kinematics fourth output.

    if nargin < 3, model=[]; end
    if ~(isnumeric(theta1) && isscalar(theta1) && ...
            isnumeric(theta2) && isscalar(theta2))
        error('jacobian_gamma_theta:InvalidInput', ...
            'theta1 and theta2 must be finite real numeric scalars in degrees.');
    end
    [~, ~, J] = coneik_kinematics([double(theta1); double(theta2)],[0;0],model);
end
