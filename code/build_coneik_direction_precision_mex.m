function output_path=build_coneik_direction_precision_mex(boost_include_dir)
%BUILD_CONEIK_DIRECTION_PRECISION_MEX Build adaptive raw-XYZ multiprecision acceleration.
% Pass a Boost header root containing boost/multiprecision/cpp_bin_float.hpp,
% or set BOOST_ROOT to that directory. Boost headers are build-time only;
% the delivered MEX needs no Boost DLL. The ordinary build_coneik_mex is intact.
    if nargin<1||isempty(boost_include_dir),boost_include_dir=getenv('BOOST_ROOT');end
    if ~isfile(fullfile(boost_include_dir,'boost','multiprecision','cpp_bin_float.hpp'))
        error('build_coneik_direction_precision_mex:BoostUnavailable',...
            'Provide the Boost header root or set BOOST_ROOT. Ordinary and symbolic solvers remain usable without this optional build.');
    end
    configurations=mex.getCompilerConfigurations('C++','Selected');
    if isempty(configurations)
        error('build_coneik_direction_precision_mex:CompilerNotConfigured','Configure a supported C++ compiler with mex -setup C++.');
    end
    source_dir=fileparts(mfilename('fullpath'));
    if ispc&&contains(configurations(1).Manufacturer,'Microsoft','IgnoreCase',true)
        flags='COMPFLAGS=$COMPFLAGS /std:c++17 /fp:strict';
    else
        flags='CXXFLAGS=$CXXFLAGS -std=c++17 -fno-fast-math -ffp-contract=off';
    end
    coneik_refresh_backends();
    mex('-R2018a','-O',flags,['-I',char(boost_include_dir)],'-outdir',source_dir,...
        '-output','coneik_direction_precision_mex',fullfile(source_dir,'native','coneik_direction_precision_mex.cpp'));
    output_path=fullfile(source_dir,['coneik_direction_precision_mex.',mexext]);
    coneik_refresh_backends();
    if nargout==0,fprintf('Built %s\n',output_path);end
end
