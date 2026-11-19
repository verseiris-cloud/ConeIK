function [q,info]=coneik_select_inverse_position(target,reference_q,model,options)
%CONEIK_SELECT_INVERSE_POSITION Select one bounded pair without motion derivatives.
% Q=CONEIK_SELECT_INVERSE_POSITION(TARGET,REFERENCE_Q,MODEL,OPTIONS) returns a
% 2-by-1 nominal joint pair in degrees. TARGET is [pitch,yaw] (degrees) or
% original [x,y,z]; XYZ is never converted to rounded direction angles.
% REFERENCE_Q defaults to [0;0]; MODEL defaults to coneik_model().
% OPTIONS: joint_bounds (2-by-2 [lo hi], default unbounded), weights (two
% nonnegative finite values, not both zero), objective ('distance', default,
% or 'ikjl_periodic'), on_error ('report', default, or 'error'), and
% direction_tolerance_deg (default 1e-8). No velocity or acceleration inputs.
% distance minimizes weighted physical angle distance. ikjl_periodic uses
% sum(w.*4.*sind(wrapped_angle_difference/2).^2), algebraically the IK.jl
% unit-circle objective; equal periodic costs prefer shorter physical travel.
% Both analytic branches are considered. Each joint's nearest admissible
% periodic representative is obtained without enumerating every turn. At an
% exact equal-cone center the free first angle is minimized analytically.
% Bounds apply strictly to returned doubles; no angle is clamped off a root.
% Large turn indices/phase loss report numerical_unresolved. The selection
% scope is this separable static objective, not trajectory or interval proof.
% Existing inverse APIs and branch order are unchanged; this is opt-in.
% Example: [q,info]=coneik_select_inverse_position([7,20],[10;-60],...
%          coneik_model(15,25),struct('objective','ikjl_periodic'));
    if nargin<2||isempty(reference_q),reference_q=[0;0];end
    if nargin<3,model=[];end
    if nargin<4,options=struct();end
    coneik_position_selection_input(target,reference_q);
    c=coneik_position_selection_options(model,options);
    target=double(target(:)');reference_q=double(reference_q(:));
    if nargout>1
            [q,info]=coneik_select_inverse_position_resolved(target,reference_q,c);
        else
            q=coneik_select_inverse_position_resolved(target,reference_q,c);
        end
end
