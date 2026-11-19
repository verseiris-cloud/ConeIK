function main(input_file, output_file, smooth_win)
%MAIN One-click batch inverse entry for this repository.
% MAIN runs the whole pipeline with defaults: read the trajectory txt,
% solve every pitch/yaw record, differentiate the recovered joint angles,
% and write the result txt next to the input.
%
%   main                       use data/sample_synthetic.dat
%   main('my_trajectory.txt')  any 3-column txt: time_s pitch_deg yaw_deg
%   main(in, out, w)           custom output path and smoothing window
%
% Editing your own file path into the EDIT HERE block at the top of the
% function body works too; command-line arguments override that block.
%
% The input txt allows one leading text header line; columns are separated
% by tabs, spaces or commas. Rows with a failed solve keep NaN angles; the
% differentiation runs per continuous valid segment exactly like
% code/batch_inverse_dat.m (unwrap in radians, central differences,
% moving-mean smoothing, default window 11 samples).
% The automatic entry uses the GPU for large batches when available.
% Small batches and systems without a usable GPU use the CPU row engine.
% The output txt has the columns
%   time_s theta1_deg theta2_deg dtheta1_deg_s dtheta2_deg_s
%   ddtheta1_deg_s2 ddtheta2_deg_s2

here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'code'));

%% ======================= EDIT HERE =====================================
%% To run your own data, change the file names below and save, then run
%% this file. Relative paths are resolved against the repository root.
%% Command-line arguments still override these values.
INPUT_FILE  = 'data/sample_synthetic.dat';          % 10001-row synthetic trajectory
OUTPUT_FILE = 'data/sample_synthetic_result.dat';   % written by the run
SMOOTH_WIN  = 11;                                   % derivative window (2 = off)
GPU_MIN_ROWS = 4096;                                % minimum rows for automatic GPU use
%% =======================================================================

if nargin<1||isempty(input_file)
    input_file=INPUT_FILE;
end
if ischar(input_file)&&~isempty(input_file)&&~isfolder(input_file)&&input_file(1)~=filesep ...
        &&~(numel(input_file)>1&&input_file(2)==':')
    cand=fullfile(here,input_file);
    if isfile(cand),input_file=cand;end
end
if nargin<2||isempty(output_file)
    output_file=OUTPUT_FILE;
end
if ischar(output_file)&&~isempty(output_file)&&output_file(1)~=filesep ...
        &&~(numel(output_file)>1&&output_file(2)==':')
    output_file=fullfile(here,output_file);
end
if nargin<3||isempty(smooth_win)
    smooth_win=SMOOTH_WIN;
end

%% read
[t,g1,g2]=read_trajectory(input_file);
n=numel(t);
fprintf('input : %s (%d records)\n',input_file,n);

%% solve (GPU for large batches when available)
gamma=[g1,g2];
model=coneik_model();
backend='cpu';
Q=[];
if n>=GPU_MIN_ROWS
    try
        % Warm up the CUDA libraries and kernel before timing the batch.
        coneik_inverse_batch_gpu(gamma(1:8,:),model);
        solve_start=tic;
        Q=coneik_inverse_batch_gpu(gamma,model);
        backend='gpu';
    catch
        % Use the CPU row engine when GPU initialization or execution fails.
    end
end
if isempty(Q)
    solve_start=tic;
    Q=coneik_inverse_batch(gamma,model);
end
solve_time=toc(solve_start);
theta1=Q(:,1);theta2=Q(:,2);
n_fail=sum(any(isnan(Q(:,1:2)),2));
fprintf('solve : %s backend, %d/%d solved, %.1f ms\n', ...
    backend,n-n_fail,n,solve_time*1e3);

%% differentiate (same scheme as code/batch_inverse_dat.m)
w=round(smooth_win);if w<2,w=1;end
dtheta1=nan(n,1);dtheta2=nan(n,1);ddtheta1=nan(n,1);ddtheta2=nan(n,1);
valid=~isnan(theta1);
seg_starts=find(diff([false;valid])==1);
seg_ends=find(diff([valid;false])==-1);
for k=1:numel(seg_starts)
    idx=seg_starts(k):seg_ends(k);
    if numel(idx)<3,continue;end
    th1=unwrap(deg2rad(theta1(idx)));th2=unwrap(deg2rad(theta2(idx)));tv=t(idx);
    if w>1
        d1=movmean(rad2deg(gradient(th1,tv)),w);d2=movmean(rad2deg(gradient(th2,tv)),w);
        dd1=movmean(gradient(d1,tv),w);dd2=movmean(gradient(d2,tv),w);
    else
        d1=rad2deg(gradient(th1,tv));d2=rad2deg(gradient(th2,tv));
        dd1=gradient(d1,tv);dd2=gradient(d2,tv);
    end
    dtheta1(idx)=d1;dtheta2(idx)=d2;ddtheta1(idx)=dd1;ddtheta2(idx)=dd2;
end

%% write
fid=fopen(output_file,'w');
if fid==-1
    error('main:CannotWrite','Cannot open the output file: %s',output_file);
end
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'time_s\ttheta1_deg\ttheta2_deg\tdtheta1_deg_s\tdtheta2_deg_s\tddtheta1_deg_s2\tddtheta2_deg_s2\n');
fprintf(fid,'%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n', ...
    [t,theta1,theta2,dtheta1,dtheta2,ddtheta1,ddtheta2].');
fprintf('output: %s\n',output_file);
end

function [t,g1,g2]=read_trajectory(filename)
%READ_TRAJECTORY 3-column time/pitch/yaw txt with an optional header.
    fid=fopen(filename,'r');
    if fid==-1
        error('main:CannotRead','Cannot open the input file: %s',filename);
    end
    cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
    % A single leading text header is allowed; a numeric first line is data.
    head=fgetl(fid);
    if ischar(head)
        vals=str2double(regexp(head,'[\s,]+','split'));
        if numel(vals)>=3&&~any(isnan(vals(1:3)))
            fseek(fid,0,'bof');
        end
    else
        fseek(fid,0,'bof');
    end
    raw=textscan(fid,'%f%f%f','Delimiter','\t ,','MultipleDelimsAsOne',true);
    t=raw{1};g1=raw{2};g2=raw{3};
    n=numel(t);
    if n==0
        error('main:NoData','No numeric rows found in: %s',filename);
    end
    if any(diff(t)<=0)
        error('main:TimeNotIncreasing','The time column must be strictly increasing.');
    end
end
