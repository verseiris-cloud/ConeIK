function fn=coneik_native_backend(kind)
% Resolve positive MEX backends once, as function handles. Missing backends
% are not cached: installing an optional MEX can be detected on the next call.
% Handles intentionally bind the initialized backend. Call the public
% coneik_refresh_backends after switching installations or changing MATLAB path.
% No mathematical result or approximate geometry is cached here.
persistent position precision direction_precision
switch kind
    case 'reset'
        position=[];precision=[];direction_precision=[];fn=[];return;
    case 'position'
        if isempty(position)&&exist('coneik_inverse_mex','file')==3,position=@coneik_inverse_mex;end
        fn=position;
    case 'precision'
        if isempty(precision)&&exist('coneik_precision_mex','file')==3,precision=@coneik_precision_mex;end
        fn=precision;
    case 'direction_precision'
        if isempty(direction_precision)&&exist('coneik_direction_precision_mex','file')==3,direction_precision=@coneik_direction_precision_mex;end
        fn=direction_precision;
    otherwise
        error('coneik_native_backend:InvalidKind','Unknown native backend kind.');
end
end
