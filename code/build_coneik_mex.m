function output_path = build_coneik_mex()
%BUILD_CONEIK_MEX Compile the optional batch position kernel using MATLAB MEX.
%   BUILD_CONEIK_MEX uses the selected C++ compiler (configure with
%   mex -setup C++), the R2018a real-double API and conservative FP flags.
%   Output is coneik_inverse_mex.<mexext> beside this function. The portable
%   native/coneik_core.hpp can also be used without MATLAB. Native status -2
%   requires the caller to run the full MATLAB precision/identity solver.
%   No compiler is downloaded or configured by this function.

    configurations = mex.getCompilerConfigurations('C++','Selected');
    if isempty(configurations)
        error('build_coneik_mex:CompilerNotConfigured', ...
            'No C++ MEX compiler is configured. Configure a supported compiler with mex -setup C++, then run build_coneik_mex again.');
    end
    source_dir = fileparts(mfilename('fullpath'));
    source_path = fullfile(source_dir,'native','coneik_mex.cpp');
    manufacturer = configurations(1).Manufacturer;
    if ispc && contains(manufacturer,'Microsoft','IgnoreCase',true)
        flags = 'COMPFLAGS=$COMPFLAGS /std:c++17 /fp:strict';
    else
        flags = 'CXXFLAGS=$CXXFLAGS -std=c++17 -fno-fast-math -ffp-contract=off';
    end
    coneik_refresh_backends();
    mex('-R2018a','-O',flags,'-outdir',source_dir,'-output','coneik_inverse_mex',source_path);
    output_path = fullfile(source_dir,['coneik_inverse_mex.',mexext]);
    coneik_refresh_backends();
    if nargout == 0
        fprintf('Built %s\n',output_path);
    end
end
