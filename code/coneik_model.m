function model = coneik_model(alpha1_deg, alpha2_deg)
%CONEIK_MODEL Construct fixed geometry for the parameterized ConeIK solver.
%   MODEL=CONEIK_MODEL(A1,A2) uses cone angles in degrees, 0<A1,A2<90.
%   MODEL=CONEIK_MODEL() selects the original 20/20 degree mechanism.
%   Geometry: r=Rx(q1)*Rz(A1)*Rx(q2)*Rz(A2)*[1;0;0].
%   Both joint angles are NOMINAL degrees. This direction-only model does
%   not describe arbitrary translational, flexible or nonconstant-axis errors.
%   Treat the returned value as immutable; rebuild it to change geometry.

    if nargin == 0
        alpha1_deg = 20; alpha2_deg = 20;
    elseif nargin ~= 2
        error('coneik_model:InvalidAngles', 'Provide both cone angles in degrees, or omit both.');
    end
    if ~valid_angle(alpha1_deg) || ~valid_angle(alpha2_deg)
        error('coneik_model:InvalidAngles', ...
            'Cone angles must be finite real numeric scalars strictly between 0 and 90 degrees.');
    end
    a = double(alpha1_deg); b = double(alpha2_deg);
    [sum_hi, sum_lo] = two_sum(a, b);
    [diff_hi, diff_lo] = two_sum(a, -b);
    sa = sind(a); sb = sind(b); ca = cosd(a); cb = cosd(b);
    % Evaluate exact binary64 sums through their two-term expansions. A
    % rounded sum of 90 must not turn a nearby geometry into a chart pole.
    cos_sum = cosd(sum_hi)*cosd(sum_lo)-sind(sum_hi)*sind(sum_lo);
    sin_sum = sind(sum_hi)*cosd(sum_lo)+cosd(sum_hi)*sind(sum_lo);
    cos_diff = cosd(diff_hi)*cosd(diff_lo)-sind(diff_hi)*sind(diff_lo);
    sin_diff = sind(diff_hi)*cosd(diff_lo)+cosd(diff_hi)*sind(diff_lo);
    model = struct('schema', 'coneik.model', 'version', 1, ...
        'alpha_deg', [a,b], 'is_default', a == 20 && b == 20, ...
        'alpha_rad', [a,b]*(pi/180), 'sin_alpha', [sa,sb], ...
        'cos_alpha', [ca,cb], 'sin_product', sa*sb, ...
        'outer_deg', sum_hi, 'inner_deg', abs(diff_hi), ...
        'sum_roundoff_deg', sum_lo, 'difference_deg', diff_hi, ...
        'difference_roundoff_deg', diff_lo, ...
        'sin_sum', sin_sum, 'cos_sum', cos_sum, ...
        'sin_difference', sin_diff, 'cos_difference', cos_diff, ...
        'half_outer_sin', sind(sum_hi/2)*cosd(sum_lo/2)+cosd(sum_hi/2)*sind(sum_lo/2), ...
        'half_inner_sin', abs(sind(diff_hi/2)*cosd(diff_lo/2)+cosd(diff_hi/2)*sind(diff_lo/2)), ...
        'key', ['coneik.model.v1:',num2hex(a),':',num2hex(b)]);
end

function yes = valid_angle(value)
    yes = isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value > 0 && value < 90;
end

function [hi, lo] = two_sum(a, b)
    hi = a+b; z = hi-a;
    lo = (a-(hi-z))+(b-z);
end
