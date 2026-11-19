function select=coneik_prepare_xyz_position_selector(model,options)
%CONEIK_PREPARE_XYZ_POSITION_SELECTOR Optional prepared XYZ-only angle selection.
% SELECT=CONEIK_PREPARE_XYZ_POSITION_SELECTOR(MODEL,OPTIONS) fixes geometry and
% distance-selection options. Q=SELECT([X,Y,Z],REFERENCE_Q) returns one 2x1
% angle pair in degrees; REFERENCE_Q defaults to [0;0]. The direction need
% not be normalized. Bounds, weights, tolerance and on_error have the same
% meaning as coneik_prepare_position_selector. Only one XYZ target and one
% output are supported. Use the existing selector for diagnostics, pitch/yaw
% or the ikjl_periodic objective.
% Optional MEX availability is sampled at construction. Recreate SELECT
% after installing/removing/replacing backends or changing MATLAB paths.
% The existing coneik_refresh_backends does not refresh this new configuration.
    if nargin<1,model=[];end
    if nargin<2,options=struct();end
    c=coneik_position_selection_options(model,options);
    if ~strcmp(c.objective,'distance')
        error('coneik_prepare_xyz_position_selector:UnsupportedObjective',...
            'This XYZ-only selector supports the distance objective.');
    end
    fast=c.tolerance>=1e-10 && ~issparse(c.bounds) && ~issparse(c.w) && ...
        ~issparse(c.model.alpha_deg) && ~issparse(c.tolerance) && ...
        exist('coneik_direction_double_mex','file')==3 && ...
        exist('coneik_position_select_mex','file')==3;
    select=@query;
    function q=query(target,reference_q)
        if nargin<2||isempty(reference_q),reference_q=[0;0];end
        coneik_position_selection_input(target,reference_q);
        if numel(target)~=3
            error('coneik_prepare_xyz_position_selector:InvalidXYZ',...
                'target must be one full real XYZ 3-vector.');
        end
        target=double(target(:)');reference_q=double(reference_q(:));
        if fast && ~issparse(reference_q)
            q=coneik_try_select_position_xyz_q(target,reference_q,c);
            if all(isfinite(q)),return;end
        end
        q=coneik_select_inverse_position_resolved(target,reference_q,c);
    end
end
