function select=coneik_prepare_position_selector(model,options)
%CONEIK_PREPARE_POSITION_SELECTOR Reuse fixed geometry and static selection options.
% SELECT=CONEIK_PREPARE_POSITION_SELECTOR(MODEL,OPTIONS) validates fixed values
% once. Q=SELECT(TARGET) or [Q,INFO]=SELECT(TARGET,REFERENCE_Q) selects a
% nominal joint pair. TARGET is pitch/yaw (degrees) or original XYZ.
% REFERENCE_Q defaults to [0;0]. Options match coneik_select_inverse_position.
% Caller-side changes do not affect SELECT. No query state or MEX handle is
% captured; coneik_refresh_backends also affects existing prepared selectors.
    if nargin<1,model=[];end
    if nargin<2,options=struct();end
    c=coneik_position_selection_options(model,options);
    select=@query;
    function [q,info]=query(target,reference_q)
        if nargin<2||isempty(reference_q),reference_q=[0;0];end
        coneik_position_selection_input(target,reference_q);
        target=double(target(:)');reference_q=double(reference_q(:));
        if nargout>1
            [q,info]=coneik_select_inverse_position_resolved(target,reference_q,c);
        else
            q=coneik_select_inverse_position_resolved(target,reference_q,c);
        end
    end
end
