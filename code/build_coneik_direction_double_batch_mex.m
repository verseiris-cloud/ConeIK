function output_path=build_coneik_direction_double_batch_mex()
%BUILD_CONEIK_DIRECTION_DOUBLE_BATCH_MEX Build optional ordinary raw-XYZ small-batch IK.
% Requires a configured C++17 MEX compiler, with strict floating point.
% No Boost or EAIK dependency. Missing MEX retains the MATLAB solver.
    configurations=mex.getCompilerConfigurations('C++','Selected');
    if isempty(configurations)
        error('build_coneik_direction_double_batch_mex:CompilerNotConfigured',...
            'Configure a supported C++ compiler with mex -setup C++.');
    end
    source_dir=fileparts(mfilename('fullpath'));
    if ispc&&contains(configurations(1).Manufacturer,'Microsoft','IgnoreCase',true)
        flags='COMPFLAGS=$COMPFLAGS /std:c++17 /fp:strict';
    else
        flags='CXXFLAGS=$CXXFLAGS -std=c++17 -fno-fast-math -ffp-contract=off';
    end
    coneik_refresh_backends();
    mex('-R2018a','-O',flags,'-outdir',source_dir,'-output','coneik_direction_double_batch_mex',...
        fullfile(source_dir,'native','coneik_direction_double_batch_mex.cpp'));
    output_path=fullfile(source_dir,['coneik_direction_double_batch_mex.',mexext]);
    coneik_refresh_backends();
    if nargout==0,fprintf('Built %s\n',output_path);end
end
