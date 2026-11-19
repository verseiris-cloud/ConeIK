function result=demo_coneik_xyz()
%DEMO_CONEIK_XYZ Direct XYZ -> two nominal angles, optional diagnostics.
% Replace example geometry, beam, bounds and reference with your own values.
% For a target point P, first use R=P-origin in the model coordinate frame.
model=coneik_model(15,25);
R=[1,.25,-.125];
reference=[10;-60];
options=struct('joint_bounds',[-180,180;-180,180],'on_error','error');

% Build once, OUTSIDE your repeated-query loop. Limits and reference choose
% one pair from the possible inverse branches. Only q is needed in the loop.
select=coneik_prepare_xyz_position_selector(model,options);
q=select(R,reference); % 2-by-1 [q1;q2], nominal degrees

% Optional diagnostic pass: does not slow calls that only use SELECT/SOLVE.
% 0.001 degree below is an ILLUSTRATIVE assumption, not sensor calibration.
assessment=coneik_assess_direction(R,q.',model,...
    struct('direction_uncertainty_deg',0.001));

% If you want BOTH candidate pairs or have independent N-by-3 targets:
solve=coneik_prepare_direction(model);
Q=solve(R); % N-by-4 [q11,q21,q12,q22]; no automatic unique selection
result=struct('model',model,'R',R,'q',q,'candidates',Q,'assessment',assessment);
fprintf('Selected nominal angles: q1 = %.12g deg, q2 = %.12g deg\n',q);
fprintf('Direction residual = %.4g deg; local inverse gain = %.4g\n',...
    assessment.direction_error_deg,assessment.inverse_gain);
fprintf('Assumed direction perturbation = %.4g deg; local joint-norm estimate = %.4g deg\n',...
    assessment.direction_uncertainty_deg,assessment.local_joint_error_estimate_deg);
end
