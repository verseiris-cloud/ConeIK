function result=demo_coneik_position_selection()
%DEMO_CONEIK_POSITION_SELECTION Beam -> one bounded pair, without derivatives.
% Replace these example cone angles and nominal bounds with your mechanism.
model=coneik_model(15,25);
reference=[10;-60];
options=struct('joint_bounds',[-180,180;-180,180],...
    'weights',[1;1],'objective','distance','on_error','error');
gamma=[7,20]; % [pitch,yaw], degrees
[q,info]=coneik_select_inverse_position(gamma,reference,model,options);

% If the sensor supplies XYZ, pass that original vector directly. This
% vector is only constructed from angles to illustrate both APIs together.
ray=[cosd(gamma(1))*cosd(gamma(2)),cosd(gamma(1))*sind(gamma(2)),-sind(gamma(1))];
[q_xyz,xyz_info]=coneik_select_inverse_position(ray,reference,model,options);

% Repeated queries: validate fixed geometry/options once. Request only Q in
% the fast loop; candidates, bounds, phase guards and FK checks still run.
select=coneik_prepare_position_selector(model,options);
q_prepared=select(ray,reference);
[q_checked,prepared_info]=select(ray,reference);
assert(isequaln(q_xyz,q_prepared)&&isequaln(q_xyz,q_checked)&&isequaln(xyz_info,prepared_info));

options.objective='ikjl_periodic';
[q_periodic,periodic_info]=coneik_select_inverse_position(ray,reference,model,options);
result=struct('model',model,'reference',reference,'gamma',gamma,'ray',ray,...
    'q',q,'info',info,'q_xyz',q_xyz,'xyz_info',xyz_info,...
    'q_periodic',q_periodic,'periodic_info',periodic_info,'q_prepared',q_prepared);
fprintf('Distance selection: q1 = %.12g deg, q2 = %.12g deg\n',q);
fprintf('XYZ direction error: %.4g deg\n',xyz_info.direction_error_deg);
end
