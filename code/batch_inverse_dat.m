function [result, stats, trajectory] = batch_inverse_dat(dat_file, smooth_win, options)
%BATCH_INVERSE_DAT 批量反解：由 .dat 方向角序列求名义摆角序列及其速度/加速度。
%
% 输入：
%   dat_file   - .dat 文件路径。每行 3 个数值：
%                第 1 列：时间 t（秒），要求严格递增；
%                第 2 列：Qv = gamma1，俯仰角（度）；
%                第 3 列：Qh = gamma2，偏航角（度）。
%                允许首行为非数值表头（如 "time Qv Qh"），自动跳过；
%                分隔符兼容制表符、空格、逗号。
%   smooth_win - 可选，多周期平滑窗长（采样点数，默认 11；1 ms 采样下即
%                跨 11 个采样周期、约 11 ms）。速度/加速度在数值微分后按该窗做
%                滑动平均（movmean），抑制量化噪声经微分放大形成的噪声带；
%                置 0 或 1 可关闭平滑。位置列（第 2、3 列）为原始反解值，不平滑。
%
% 输出：
%   result - N×7 矩阵，每行与输入数据行一一对应：
%            [t, theta1, theta2, dtheta1/dt, dtheta2/dt, d2theta1/dt2, d2theta2/dt2]
%            单位：秒、度、度、度/秒、度/秒、度/秒^2、度/秒^2。
%            theta1、theta2 为名义摆角（inverse_motion_analysis 第一组逆解分支，
%            即经 theta2 = theta2^r - theta1^r 变换后的名义量，与反解 GUI 显示一致）。
%            速度为 theta 序列对时间的中心差分（gradient）再经多周期滑动平均；
%            加速度为对平滑后速度再做一次同样的差分与滑动平均；再次滤波后，
%            它不是输出速度列的严格未滤波导数。微分前在弧度下 unwrap，避免 ±180° 角度
%            回绕污染导数。超出可达域或反解失败的行：各列均为 NaN。
%   stats  - 统计结构体：n_total（总行数）、n_ok（成功行数）、
%            n_fail（失败行数）、fail_rows（失败行号向量）。
%
% 模式：无 options 的旧调用保持历史解析、固定默认模型及差分定义，
% 位置计算默认自动使用 C++（可用时），浮点末位可能与 MATLAB 有差异。
% 显式 options 必须带 mode='position'/'derivative'/'continuous'；三种模式
% 都可用 model=coneik_model(alpha1,alpha2)，锥角示例不是固定限制。
% position：独立逐点第一分支，N×7 中第4:7列为NaN，trajectory=[]。
% derivative：独立位置 + 分段 unwrap/gradient/movmean，非可执行轨迹保证。
% continuous：q/v/a 来自同一曲线并接受约束检查，详见 coneik_inverse_trajectory。
% position/continuous 的 smooth_win 只接受 []/0/1；derivative 默认11。
% 显式模式严格解析 DAT，NaN/不可达方向保留为缺口。position/derivative
% 只接受 mode、model、position_backend（默认auto）；连续约束不能静默忽略。
% stats.n_ok/n_fail 在独立模式按位置成败统计；derivative_available 另报导数。
% 常规 ASCII 三列 DAT 按块批量解析；特殊格式回退历史逐行解析，保留其语义。

    trajectory=[];
    if nargin<2,smooth_win=[];end
    if nargin>=3 && ~isempty(options)
        if ~(isstruct(options) && isscalar(options))
            error('batch_inverse_dat:InvalidOptions','options must be a scalar struct.');
        end
        if ~isempty(fieldnames(options))
            if any(isfield(options,{'trajectory_selection_policy','trajectory_selection_weights'}))
                error('coneik_inverse_trajectory:RemovedFeature',...
                    'Trajectory smoothness ranking is temporarily removed. Omit trajectory_selection_policy and trajectory_selection_weights.');
            end
            if ~isfield(options,'mode')
                error('batch_inverse_dat:ModeRequired','Explicit options require mode=position, derivative or continuous.');
            end
            mode=options.mode;
            if isstring(mode)&&isscalar(mode)&&~ismissing(mode),mode=char(mode);end
            if ~(ischar(mode)&&isrow(mode)&&any(strcmp(mode,{'position','derivative','continuous'})))
                error('batch_inverse_dat:InvalidOptions','mode must be position, derivative or continuous.');
            end
            if ~strcmp(mode,'derivative') && ...
                    ~(isempty(smooth_win)||(isnumeric(smooth_win)&&isreal(smooth_win)&&...
                    isscalar(smooth_win)&&isfinite(smooth_win)&&any(smooth_win==[0,1])))
                error('batch_inverse_dat:IncompatibleSmoothing','Position/continuous mode requires smooth_win=[], 0 or 1.');
            end
            if strcmp(mode,'continuous')
                if ~isfield(options,'position_backend'),options.position_backend='auto';end
                [times,angles,line_numbers]=read_continuous_dat(dat_file);
                [result,stats,trajectory]=coneik_inverse_trajectory(times,angles,options);
                stats.input_line_numbers=line_numbers;
            else
                if any(~ismember(fieldnames(options),{'mode','model','position_backend'}))
                    error('batch_inverse_dat:InvalidOptions','Independent modes only accept mode, model and position_backend; use continuous for trajectory constraints.');
                end
                if strcmp(mode,'derivative')
                    if isempty(smooth_win),smooth_win=11;end
                    if ~(isnumeric(smooth_win)&&isreal(smooth_win)&&isscalar(smooth_win)&&...
                            isfinite(smooth_win)&&smooth_win>=0&&smooth_win<=flintmax/2)
                        error('batch_inverse_dat:InvalidSmoothing','Derivative smoothing window must be a finite nonnegative scalar.');
                    end
                end
                [times,angles,line_numbers]=read_continuous_dat(dat_file);
                model=[];if isfield(options,'model'),model=options.model;end
                backend='auto';if isfield(options,'position_backend'),backend=options.position_backend;end
                [Q,info]=coneik_inverse_batch(angles,model,struct('backend',backend,'on_error','report'));
                n=numel(times);result=[times,Q(:,1:2),nan(n,4)];
                if strcmp(mode,'derivative')
                    result(:,4:7)=independent_derivatives(times,Q(:,1:2),smooth_win);
                end
                available=all(isfinite(result(:,4:7)),2);
                stats=struct('n_total',n,'n_ok',nnz(info.success),'n_fail',nnz(~info.success),...
                    'fail_rows',find(~info.success)','mode',mode,'model_key',info.model_key,...
                    'input_line_numbers',line_numbers,'position_info',info,...
                    'derivatives_computed',strcmp(mode,'derivative'),'derivative_available',available,...
                    'derivative_unavailable_rows',find(~available)',...
                    'continuous_motion_verified',false,'branch','first independent analytic branch',...
                    'scope','independent positions; derivatives are sampled analysis only; no continuous executability guarantee');
                if strcmp(mode,'derivative'),stats.smoothing_window=max(1,round(double(smooth_win)));end
            end
            return;
        end
    end

    % ---------- 读取与解析 ----------
    fid = fopen(dat_file, 'r');
    if fid == -1
        error('batch_inverse_dat: 无法打开文件：%s', dat_file);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    [t,g,~,parsed,resume]=coneik_read_dat_ascii(fid,false);
    if parsed
        g1=g(:,1);g2=g(:,2);
    else
        if resume.available && fseek(fid,resume.offset,'bof')==0
            g1=g(:,1);g2=g(:,2);
        else
            fseek(fid,0,'bof');
            t = [];  g1 = [];  g2 = [];
        end
    end
    while ~parsed
        ln = fgetl(fid);
        if ~ischar(ln)
            break;
        end
        ln = strtrim(ln);
        if isempty(ln)
            continue;
        end
        vals = str2double(regexp(ln, '[\s,]+', 'split'));   % 兼容制表符/空格/逗号
        if numel(vals) < 3 || any(isnan(vals(1:3)))
            continue;                                       % 跳过表头及非数值行
        end
        t(end+1, 1)  = vals(1);                             %#ok<AGROW>
        g1(end+1, 1) = vals(2);                             %#ok<AGROW>
        g2(end+1, 1) = vals(3);                             %#ok<AGROW>
    end

    n = numel(t);
    if n == 0
        error('batch_inverse_dat: 文件中没有有效数据行：%s', dat_file);
    end
    if any(diff(t) <= 0)
        error('batch_inverse_dat: 时间列必须严格递增。');
    end

    % ---------- 逐行反解名义摆角（第一组分支） ----------
    [positions,position_info]=coneik_inverse_batch([g1,g2],[],struct('backend','auto'));
    theta1=positions(:,1);theta2=positions(:,2);
    fail_rows=find(~position_info.success)';

    % ---------- 对 theta 序列做时间数值微分（逐有效连续段） ----------
    % 速度 = movmean(gradient(theta), w)；加速度 = movmean(gradient(平滑后速度), w)
    if nargin < 2 || isempty(smooth_win)
        smooth_win = 11;               % 默认多周期平滑窗长（采样点数）
    end
    w = round(smooth_win);
    if w < 2
        w = 1;                         % 0/1 表示关闭平滑
    end

    dtheta1  = nan(n, 1);  dtheta2  = nan(n, 1);
    ddtheta1 = nan(n, 1);  ddtheta2 = nan(n, 1);

    valid = ~isnan(theta1);
    seg_starts = find(diff([false; valid]) == 1);
    seg_ends   = find(diff([valid; false]) == -1);
    for k = 1:numel(seg_starts)
        idx = seg_starts(k):seg_ends(k);
        if numel(idx) < 3
            continue;                % 段太短无法中心差分，导数保持 NaN
        end
        th1 = unwrap(deg2rad(theta1(idx)));       % 弧度下 unwrap 防 ±180° 跳变
        th2 = unwrap(deg2rad(theta2(idx)));
        tv  = t(idx);
        if w > 1
            d1  = movmean(rad2deg(gradient(th1, tv)), w);
            d2  = movmean(rad2deg(gradient(th2, tv)), w);
            dd1 = movmean(gradient(d1, tv), w);
            dd2 = movmean(gradient(d2, tv), w);
        else
            d1  = rad2deg(gradient(th1, tv));
            d2  = rad2deg(gradient(th2, tv));
            dd1 = gradient(d1, tv);
            dd2 = gradient(d2, tv);
        end
        dtheta1(idx)  = d1;   dtheta2(idx)  = d2;
        ddtheta1(idx) = dd1;  ddtheta2(idx) = dd2;
    end

    result = [t, theta1, theta2, dtheta1, dtheta2, ddtheta1, ddtheta2];
    stats = struct('n_total', n, ...
                   'n_ok', n - numel(fail_rows), ...
                   'n_fail', numel(fail_rows), ...
                   'fail_rows', fail_rows);
    stats.position_backend=position_info.backend_used;
    stats.native_fast_count=position_info.native_fast_count;
    stats.matlab_row_count=position_info.matlab_row_count;
end

function [t,g,lines]=read_continuous_dat(filename)
% A single initial textual header is allowed. No data row is silently lost.
    fid=fopen(filename,'r');
    if fid<0, error('batch_inverse_dat:OpenFailed','Cannot open input file.'); end
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    [t,g,lines,parsed,resume]=coneik_read_dat_ascii(fid,true);
    line_number=0;first=true;
    if ~parsed
        if resume.available && fseek(fid,resume.offset,'bof')==0
            line_number=resume.line_count;first=resume.first;
        else
            fseek(fid,0,'bof');t=zeros(0,1);g=zeros(0,2);lines=zeros(0,1);
        end
    end
    while ~parsed
        line=fgetl(fid); if ~ischar(line), break; end
        line_number=line_number+1; line=strtrim(line);
        if isempty(line), continue; end
        if ~isempty(regexp(line,'(^,)|(,$)|(,\s*,)','once'))
            error('batch_inverse_dat:MalformedRecord','Empty field at line %d.',line_number);
        end
        tokens=regexp(line,'[\s,]+','split'); values=str2double(tokens);
        numeric=~isnan(values) | strcmpi(tokens,'nan');
        if first && numel(tokens)==3 && all(~numeric)
            first=false; continue;
        end
        first=false;
        if numel(values)~=3 || any(~numeric)
            error('batch_inverse_dat:MalformedRecord',...
                'Malformed input record at line %d; expected time,pitch,yaw.',line_number);
        end
        if ~isfinite(values(1))
            error('batch_inverse_dat:InvalidTime','Nonfinite time at line %d.',line_number);
        end
        t(end+1,1)=values(1); g(end+1,:)=values(2:3); lines(end+1,1)=line_number; %#ok<AGROW>
    end
    if isempty(t), error('batch_inverse_dat:EmptyInput','No numeric records in input file.'); end
    if any(diff(t)<=0)
        error('batch_inverse_dat:InvalidTime','Time must be strictly increasing.');
    end
end

function derivatives=independent_derivatives(t,q,smooth_win)
% Same segmented unwrap/gradient/movmean definition as the legacy branch.
    n=numel(t);derivatives=nan(n,4);w=max(1,round(double(smooth_win)));
    valid=all(isfinite(q),2);
    starts=find(diff([false;valid])==1);ends=find(diff([valid;false])==-1);
    for k=1:numel(starts)
        idx=starts(k):ends(k);if numel(idx)<3,continue;end
        for axis=1:2
            theta=unwrap(deg2rad(q(idx,axis)));
            v=rad2deg(gradient(theta,t(idx)));
            if w>1,v=movmean(v,w);end
            a=gradient(v,t(idx));if w>1,a=movmean(a,w);end
            derivatives(idx,axis)=v;derivatives(idx,axis+2)=a;
        end
    end
end
