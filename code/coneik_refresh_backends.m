function coneik_refresh_backends()
%CONEIK_REFRESH_BACKENDS Forget initialized optional native backend handles.
% Call after installing/replacing/removing a MEX, or switching installations
% on the MATLAB path. The next native solve resolves the backend again.
% This also releases the precision and ordinary-XYZ MEX geometry caches.
% Existing mathematical results and MATLAB model/state caches are not modified.
% backend='matlab' never resolves or calls a native backend.
coneik_native_backend('reset');
clear coneik_inverse_mex coneik_precision_mex coneik_direction_precision_mex coneik_direction_double_mex coneik_direction_double_batch_mex coneik_direction_resolved
rehash;
end
