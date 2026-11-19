function [result, gamma1, gamma2] = motion_analysis(theta1, theta2, model)
% Optional final MODEL is constructed by coneik_model; omitted means 20/20.
%MOTION_ANALYSIS Parameterized forward kinematics (default 20/20 degrees), nominal angles.
%   Inputs are finite real numeric scalars in degrees. Outputs preserve the
%   existing interface: unit direction, principal pitch, principal yaw.
%   The shared analytic chain already includes the nominal-angle convention.

    if nargin < 3, model=[]; end
    if ~(isnumeric(theta1) && isscalar(theta1) && ...
            isnumeric(theta2) && isscalar(theta2))
        error('motion_analysis:InvalidInput', ...
            'theta1 and theta2 must be finite real numeric scalars in degrees.');
    end
    if nargout < 2
        result = coneik_kinematics([double(theta1); double(theta2)],[0;0],model);
    else
        [result, gamma] = coneik_kinematics([double(theta1); double(theta2)],[0;0],model);
        gamma1 = gamma(1); gamma2 = gamma(2);
    end
end
