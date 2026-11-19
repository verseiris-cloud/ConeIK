function build_coneik_position_select_mex()
%BUILD_CONEIK_POSITION_SELECT_MEX Optional ordinary XYZ distance post-selector.
root=fileparts(mfilename('fullpath'));source=fullfile(root,'native','coneik_position_select_mex.cpp');
config=mex.getCompilerConfigurations('C++','Selected');
if isempty(config),error('coneik_position_select_mex:Compiler','Configure a C++17 compiler with mex -setup C++.');end
if ispc&&contains(config.Manufacturer,'Microsoft')
    flags={'COMPFLAGS=$COMPFLAGS /std:c++17 /O2 /fp:strict'};
else
    flags={'CXXFLAGS=$CXXFLAGS -std=c++17 -O3 -fno-fast-math -ffp-contract=off'};
end
mex('-R2018a',flags{:},source,'-outdir',root,'-output','coneik_position_select_mex');
coneik_refresh_backends;
end
