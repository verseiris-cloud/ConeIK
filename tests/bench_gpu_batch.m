%BENCH_GPU_BATCH CPU versus GPU timing for the ConeIK batch core.
% Run from anywhere:
%   addpath('tests'); bench_gpu_batch
% Requires the Parallel Computing Toolbox. Sections:
%   1. ordinary pool (524,288 seeded random interior directions)
%      - correctness: GPU vs CPU candidate agreement
%      - batch-size sweep 8..524,288, median wall time, speedup
%      - CPU single-target cost
%   2. fallback rows (center, boundary, unreachable, NaN, zero)
%   3. measured command sequence, only when the environment variable
%      CONEIK_H1_DATA points to a 3-column time/pitch/yaw .dat file that
%      is kept outside this repository
% Results are printed and written to tests/results/bench_gpu_batch.txt.

root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code'));
outdir=fullfile(root,'tests','results');
if ~isfolder(outdir),mkdir(outdir);end
logfile=fullfile(outdir,'bench_gpu_batch.txt');
if isfile(logfile),delete(logfile);end
diary(logfile);

try, parallel.gpu.enableCUDAForwardCompatibility(true); catch, end
g=gpuDevice;
fprintf('bench_gpu_batch  %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf('GPU: %s, CC %s, %.1f GB free of %.1f GB\n',g.Name,g.ComputeCapability,g.AvailableMemory/2^30,g.TotalMemory/2^30);

model=coneik_model(20,20);
solve_cpu=coneik_prepare_direction(model);
solve_gpu=coneik_prepare_direction_gpu(model);

%% 1. ordinary pool
rng(0);
target=524288;R=zeros(0,3);
while size(R,1)<target
    v=randn(500000,3);
    beta=atan2(hypot(v(:,2),v(:,3)),v(:,1))*180/pi;
    keep=beta>1&beta<39;              % 1 deg margin from both boundaries
    R=[R;v(keep,:)];
end
R=R(1:target,:);

Qc=solve_cpu(R);Qg=solve_gpu(R);
d=mod(Qg-Qc+180,360)-180;
fprintf('\n--- ordinary pool: %d rows ---\n',target);
fprintf('GPU rows finite      : %d / %d\n',sum(all(isfinite(Qg),2)),target);
fprintf('max |Qgpu - Qcpu|    : %.3e deg (wrapped)\n',max(abs(d(:))));

sizes=[8 32 128 512 1024 4096 8192 32768 131072 524288];
reps =[30 30 20 10 10 5 5 5 5 3];
fprintf('\n%6s %10s %12s %10s %12s %9s\n','N','cpu us','cpu us/tgt','gpu us','gpu us/tgt','speedup');
for i=1:numel(sizes)
    B=sizes(i);Rb=R(1:B,:);
    solve_gpu(Rb);wait(g);            % warm-up includes kernel build
    tc=zeros(1,reps(i));
    for k=1:reps(i),tic;solve_cpu(Rb);tc(k)=toc;end
    tg=zeros(1,reps(i));
    for k=1:reps(i),tic;solve_gpu(Rb);wait(g);tg(k)=toc;end
    mc=median(tc)*1e6;mg=median(tg)*1e6;
    fprintf('%6d %10.1f %12.3f %10.1f %12.3f %8.2fx\n',B,mc,mc/B,mg,mg/B,mc/mg);
end

Rs=R(1:200,:);ts=zeros(1,200);
for k=1:200,tic;solve_cpu(Rs(k,:));ts(k)=toc;end
fprintf('cpu single-target    : median %.2f us, p95 %.2f us\n',median(ts)*1e6,prctile(ts,95)*1e6);

%% 2. fallback rows
fprintf('\n--- fallback rows (GPU decline -> CPU row engine) ---\n');
Rspec=[1 0 0;                       % exact center, equal cones
       cosd(39.5) sind(39.5) 0;     % interior near the outer boundary
       cosd(40) sind(40) 0;         % effectively on the outer boundary
       cosd(50) sind(50) 0;         % outside the reachable annulus
       NaN NaN NaN;
       0 0 0];
Qs=solve_gpu(Rspec);
labels={'center','interior 39.5 deg','boundary 40 deg','unreachable 50 deg','NaN row','zero row'};
for k=1:size(Rspec,1)
    if all(isfinite(Qs(k,:)))
        fprintf('%-20s -> %10.4f %10.4f %10.4f %10.4f\n',labels{k},Qs(k,:));
    else
        fprintf('%-20s -> NaN (declined or rejected, no throw)\n',labels{k});
    end
end
assert(max(abs(Qs(1,:)-[0 180 0 -180]))<1e-9,'center row must recover the identity solution');

%% 3. measured command sequence (optional)
h1path=getenv('CONEIK_H1_DATA');
if ~isempty(h1path)&&isfile(h1path)
    seq=readmatrix(h1path,'NumHeaderLines',1);
    gammaIn=seq(:,[2 3]);
    p=gammaIn(:,1);y=gammaIn(:,2);cp=cosd(p);
    R1=[cp.*cosd(y),cp.*sind(y),-sind(p)];
    fprintf('\n--- measured sequence: %d records ---\n',size(R1,1));
    Qc1=solve_cpu(R1);Qg1=solve_gpu(R1);
    d1=mod(Qg1-Qc1+180,360)-180;
    fprintf('GPU rows finite      : %d / %d\n',sum(all(isfinite(Qg1),2)),size(R1,1));
    fprintf('max |Qgpu - Qcpu|    : %.3e deg (wrapped)\n',max(abs(d1(:))));
    solve_gpu(R1);wait(g);
    tc1=[];tic;Qc1=solve_cpu(R1);tc1=toc;
    tg1=[];tic;Qg1=solve_gpu(R1);wait(g);tg1=toc;
    fprintf('cpu vectorized batch : %8.2f ms (%.3f us/tgt)\n',tc1*1e3,tc1*1e6/size(R1,1));
    fprintf('gpu batch            : %8.2f ms (%.3f us/tgt), speedup %.2fx\n',tg1*1e3,tg1*1e6/size(R1,1),tc1/tg1);
    tic;Qrow=coneik_inverse_batch(gammaIn,model);trow=toc;
    fprintf('cpu row engine       : %8.2f ms (%.3f us/tgt)\n',trow*1e3,trow*1e6/size(R1,1));
    drow=mod(Qrow-Qc1+180,360)-180;
    fprintf('row engine vs cpu    : max %.3e deg\n',max(abs(drow(:))));
else
    fprintf('\n(measured sequence skipped: CONEIK_H1_DATA not set)\n');
end

fprintf('\nlog written to %s\n',logfile);
diary off;
