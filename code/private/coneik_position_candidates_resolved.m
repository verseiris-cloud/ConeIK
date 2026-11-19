function [Q,meta]=coneik_position_candidates_resolved(gamma,options,cache,model,sort_lifts)
% Private candidate work; MODEL and GAMMA are validated by the caller.
% Only trajectory search may omit local ordering: it sorts the complete
% layer by cost, branch, physical angles and parent before retaining paths.
    if nargin<5,sort_lifts=true;end
    model_key='';if ~isempty(model),model_key=model.key;end
    previous=[];
    if isfield(options,'previous_q') && ~isempty(options.previous_q)
        previous=vector2(options.previous_q,'previous_q');
    end
    bounds=[-Inf,Inf;-Inf,Inf];
    if isfield(options,'joint_bounds')
        bounds=options.joint_bounds;
        if ~(isnumeric(bounds)&&isreal(bounds)&&isequal(size(bounds),[2,2])&&...
                ~any(isnan(bounds(:)))&&all(bounds(:,1)<=bounds(:,2))&&...
                ~any(bounds(:,1)==Inf)&&~any(bounds(:,2)==-Inf))
            error('coneik_position_candidates:InvalidOptions','joint_bounds must have valid [lower upper] rows.');
        end
        bounds=double(bounds);
    end
    max_candidates=integer_option(options,'max_candidates',256,1);
    winding_radius=integer_option(options,'winding_radius',1,0);
    center_q1=[];
    if isfield(options,'center_q1') && ~isempty(options.center_q1)
        center_q1=options.center_q1;
        if ~(isnumeric(center_q1)&&isreal(center_q1)&&isvector(center_q1)&&all(isfinite(center_q1(:))))
            error('coneik_position_candidates:InvalidOptions','center_q1 must contain finite real seed angles.');
        end
        center_q1=double(center_q1(:)');
    end
    if isempty(cache)
        if isempty(model)||model.is_default
            [q11,q21,q12,q22,position_info]=inverse_motion_analysis(gamma(1),gamma(2));
        else
            [q11,q21,q12,q22,position_info]=coneik_inverse_position_resolved(gamma(1),gamma(2),model);
        end
        base=[q11,q12;q21,q22];
        cache=struct('input_gamma',gamma,'base_candidates',base,'position_info',position_info);
        if isfield(options,'model'), cache.model_key=model_key; end
    else
        if ~(isstruct(cache)&&isscalar(cache)&&all(isfield(cache,{'input_gamma','base_candidates','position_info'}))&&...
                isequal(cache.input_gamma,gamma)&&isnumeric(cache.base_candidates)&&isreal(cache.base_candidates)&&...
                isequal(size(cache.base_candidates),[2,2])&&all(isfinite(cache.base_candidates(:)))&&...
                isstruct(cache.position_info)&&isscalar(cache.position_info)&&isfield(cache.position_info,'status')&&...
                any(strcmp(cache.position_info.status,{'regular','center','inner_boundary','outer_boundary'})))
            error('coneik_position_candidates:InvalidCache','Cache must come from this helper for the identical gamma.');
        end
        expected=model;if isempty(expected),expected=coneik_resolve_model([]);end
        if isfield(cache,'model_key')
            matches=isequal(cache.model_key,expected.key);
        else
            matches=expected.is_default;
        end
        if ~matches
            error('coneik_position_candidates:InvalidCache','Cache was created for a different model.');
        end
        base=cache.base_candidates; position_info=cache.position_info;
    end
    position_tolerance=1e-9;
    is_center=strcmp(position_info.status,'center');
    if is_center
        default_q1=min(max(0,bounds(1,1)),bounds(1,2));
        seeds=[default_q1,center_q1];
        if ~isempty(previous), seeds=[previous(1),seeds]; end
        seeds=rem(seeds,360); seeds(seeds>=180)=seeds(seeds>=180)-360;
        seeds(seeds< -180)=seeds(seeds< -180)+360;
        seeds=unique(seeds,'stable');
        positions=[seeds;180*ones(size(seeds))];
        branch_ids=zeros(1,numel(seeds));
    elseif any(strcmp(position_info.status,{'inner_boundary','outer_boundary'}))
        positions=base(:,1); branch_ids=1;
    else
        positions=base; branch_ids=[1,2];
    end
    Q=zeros(2,0); branches=zeros(1,0);
    for k=1:size(positions,2)
        reference=positions(:,k);
        if ~isempty(previous), reference=previous; end
        x=lifts(positions(1,k),bounds(1,:),reference(1),winding_radius,max_candidates,position_tolerance,sort_lifts);
        y=lifts(positions(2,k),bounds(2,:),reference(2),winding_radius,max_candidates,position_tolerance,sort_lifts);
        count=numel(x)*numel(y);
        if count>max_candidates-size(Q,2)
            error('coneik_position_candidates:CandidateLimitExceeded',...
                'Candidate set exceeds max_candidates=%d. Raise the budget or narrow joint_bounds; no candidates were silently truncated.',max_candidates);
        end
        for i=1:numel(x)
            for j=1:numel(y)
                candidate=[x(i);y(j)];
                if is_center && any(all(Q==candidate,1)), continue; end
                Q(:,end+1)=candidate; branches(end+1)=branch_ids(k); %#ok<AGROW>
            end
        end
    end
    meta=struct('branch',branches,'position_info',position_info,...
        'candidate_count',size(Q,2),'position_enumeration_complete',~is_center&&all(isfinite(bounds(:))),...
        'center_family_sampled',is_center,'periodic_lifts_limited',any(~isfinite(bounds(:))),...
        'winding_radius',winding_radius,'max_candidates',max_candidates,...
        'joint_bounds',bounds,'position_tolerance',position_tolerance,...
        'coordinate_system','nominal degrees','base_candidates',base,'cache',cache);
    if isfield(options,'model'), meta.model_key=model_key; end
end

function values=lifts(base,bounds,reference,radius,budget,tol,sort_lifts)
    low=-Inf; high=Inf;
    if isfinite(bounds(1)), low=ceil((bounds(1)-base-tol)/360); end
    if isfinite(bounds(2)), high=floor((bounds(2)-base+tol)/360); end
    if low>high, values=zeros(1,0); return; end
    nearest=round((reference-base)/360);
    nearest=min(max(nearest,low),high);
    if ~all(isfinite(bounds))
        low=max(low,nearest-radius); high=min(high,nearest+radius);
    end
    count=high-low+1;
    if ~isfinite(count) || count>budget
        error('coneik_position_candidates:CandidateLimitExceeded',...
            'A periodic axis exceeds max_candidates=%d; no silent truncation is performed.',budget);
    end
    if ~all(isfinite([low,high])) || any(abs([low,high])>flintmax/2)
        error('coneik_position_candidates:NumericalUnresolved','Periodic turn indices exceed reliable integer resolution.');
    end
    turns=low:high;
    values=base+360*turns;
    values=values(values>=bounds(1)-tol & values<=bounds(2)+tol);
    if sort_lifts
        [~,order]=sortrows([abs(values(:)-reference),values(:)],[1,2]);
        values=values(order);
    end
end

function value=integer_option(options,name,default,minimum)
    value=default;
    if isfield(options,name), value=options.(name); end
    if ~(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=minimum&&value==fix(value)&&value<=flintmax/2)
        error('coneik_position_candidates:InvalidOptions','%s must be an integer >= %d.',name,minimum);
    end
    value=double(value);
end

function value=vector2(value,name)
    if ~(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==2&&all(isfinite(value(:))))
        error('coneik_position_candidates:InvalidInput','%s must be a finite real 2-vector.',name);
    end
    value=double(value(:));
end
