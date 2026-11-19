function [Q,info]=coneik_direction_resolved(R,model,policy,tolerance,compact)
% Internal XYZ engine with validated numeric input, immutable geometry/options.
    persistent ordinary_mex
    if nargin<5,compact=false;end
    % Fast success is confined to comfortably interior single targets.
    % All failures, boundary cases and tight tolerances retain the old path.
    if size(R,1)==1 && tolerance>=1e-10
        accepted=false;route='vector_half_angle';
        % Exact structural cases cannot enter the ordinary native path.
        identity=(R(2)==0 && R(3)==0) || ...
            (R(1)==0 && model.outer_deg==90 && model.sum_roundoff_deg==0);
        if ~identity
            % Cache only positive availability, never a miss or a target.
            % Public refresh clears this function's persistent state.
            if isempty(ordinary_mex) && exist('coneik_direction_double_mex','file')==3
                ordinary_mex=true;
            end
            if ~isempty(ordinary_mex)
                if compact || nargout<2
                    q=coneik_direction_double_mex(R,model.alpha_deg,tolerance);
                else
                    [q,residual]=coneik_direction_double_mex(R,model.alpha_deg,tolerance);
                end
                accepted=all(isfinite(q)); % Decline is all NaN; MEX validates both branches.
                if accepted,route='vector_cpp_double';end
            end
        end
        if ~accepted,[q,residual,accepted]=ordinary_single(R,model,tolerance);end
        if accepted
            Q=q;
            if compact
                info=true; % Internal ordinary-success token; no public diagnostics requested.
            elseif nargout>1
                info=struct('model_key',model.key,'input_representation','direction_xyz',...
                    'success',true,'solver_success',true,'status',{{'regular'}},...
                    'route',{{route}},'precision_dps',16,'used_precision',false,...
                    'branch_count',2,'direction_error_deg',residual,'error_id',{{''}},'error_message',{{''}});
            end
            return;
        end
    end
    if nargout<2
        Q=values_only(R,model,policy,tolerance);
        return;
    end
    n=size(R,1);Q=nan(n,4);
    info=struct('model_key',model.key,'input_representation','direction_xyz',...
        'success',false(n,1),'solver_success',false(n,1),'status',{repmat({'numerical_unresolved'},n,1)},...
        'route',{repmat({''},n,1)},'precision_dps',nan(n,1),'used_precision',false(n,1),...
        'branch_count',zeros(n,1),'direction_error_deg',nan(n,2),...
        'error_id',{repmat({''},n,1)},'error_message',{repmat({''},n,1)});
    % Share model constants and array arithmetic across ordinary rows. The
    % strict/tiny-tolerance and difficult paths below keep the original code.
    fast=false(n,1);
    if n>=8 && tolerance>=1e-10
        if n<=16
            [fastQ,fast,fastResidual,native]=small_batch(R,model,tolerance);
        else
            [fastQ,fast,fastResidual]=ordinary_batch(R,model,tolerance);
        end
        Q(fast,:)=fastQ(fast,:);info.success(fast)=true;info.solver_success(fast)=true;
        info.status(fast)={'regular'};info.route(fast)={'vector_half_angle'};
        if n<=16,info.route(native)={'vector_cpp_double_batch'};end
        info.precision_dps(fast)=16;info.branch_count(fast)=2;
        info.direction_error_deg(fast,:)=fastResidual(fast,:);
    end
    for k=1:n
        if fast(k),continue;end
        try
            original=R(k,:);scale=max(abs(original));
            if ~all(isfinite(original)) || scale==0
                error('coneik_inverse_direction:InvalidInput','Each direction must be finite and nonzero.');
            end
            scaled=original/scale;r=scaled/sqrt(sum(scaled.^2));
            [q,status,dps,route]=solve(original,r,model);
            info.solver_success(k)=true;info.route{k}=route;
            info.precision_dps(k)=dps;info.used_precision(k)=dps>16;
            info.branch_count(k)=2;
            if strcmp(status,'center'),info.branch_count(k)=Inf;
            elseif ~strcmp(status,'regular'),info.branch_count(k)=1;end
            residual=assess(q,r,model);
            info.direction_error_deg(k,:)=residual;
            if any(~isfinite(residual)) || any(residual>tolerance)
                error('coneik_inverse_direction:ResidualRejected','Returned angles exceed the direction tolerance.');
            end
            Q(k,:)=q;info.success(k)=true;info.status{k}=status;
        catch ex
            switch ex.identifier
                case 'coneik_inverse_direction:InvalidInput',status='invalid_input';
                case 'inverse_motion_analysis:Unreachable',status='unreachable';
                case 'inverse_motion_analysis:PrecisionUnavailable',status='precision_unavailable';
                case 'inverse_motion_analysis:NumericalUnresolved',status='numerical_unresolved';
                case 'coneik_inverse_direction:ResidualRejected',status='residual_rejected';
                otherwise,rethrow(ex);
            end
            info.status{k}=status;info.error_id{k}=ex.identifier;info.error_message{k}=ex.message;
            if strcmp(policy,'error'),error(ex.identifier,'Direction row %d: %s',k,ex.message);end
        end
    end
end

function Q=values_only(R,model,policy,tolerance)
% Omit unused diagnostic storage, never geometric or residual checks.
% Keep the same fast-row threshold and original-double fallback sequence.
    n=size(R,1);Q=nan(n,4);fast=false(n,1);
    if n>=8 && tolerance>=1e-10
        if n<=16
            [Q,fast]=small_batch(R,model,tolerance);
        else
            [Q,fast]=ordinary_batch(R,model,tolerance);
        end
    end
    for k=1:n
        if fast(k),continue;end
        try
            original=R(k,:);scale=max(abs(original));
            if ~all(isfinite(original)) || scale==0
                error('coneik_inverse_direction:InvalidInput','Each direction must be finite and nonzero.');
            end
            scaled=original/scale;r=scaled/sqrt(sum(scaled.^2));
            q=solve(original,r,model);
            residual=assess(q,r,model);
            if any(~isfinite(residual)) || any(residual>tolerance)
                error('coneik_inverse_direction:ResidualRejected','Returned angles exceed the direction tolerance.');
            end
            Q(k,:)=q;
        catch ex
            switch ex.identifier
                case {'coneik_inverse_direction:InvalidInput','inverse_motion_analysis:Unreachable',...
                        'inverse_motion_analysis:PrecisionUnavailable','inverse_motion_analysis:NumericalUnresolved',...
                        'coneik_inverse_direction:ResidualRejected'}
                    if strcmp(policy,'error'),error(ex.identifier,'Direction row %d: %s',k,ex.message);end
                otherwise,rethrow(ex);
            end
        end
    end
end

function [Q,fast,residual,native]=small_batch(R,m,tolerance)
% Caller confines this optional path to N=8..16 and ordinary tolerance.
% Positive availability only; public refresh clears this local cache.
    persistent batch_mex
    native=false(size(R,1),1);
    if isempty(batch_mex) && exist('coneik_direction_double_batch_mex','file')==3
        batch_mex=true;
    end
    if isempty(batch_mex)
        [Q,fast,residual]=ordinary_batch(R,m,tolerance);return;
    end
    if nargout>2
        [Q,native,residual]=coneik_direction_double_batch_mex(R,m.alpha_deg,tolerance);
    else
        [Q,native]=coneik_direction_double_batch_mex(R,m.alpha_deg,tolerance);
    end
    fast=native;if all(native),return;end
    % Retry every declined row with the original MATLAB array solver before
    % the established row-wise precision/error path. Original XYZ is retained.
    rows=find(~native);
    [fallback,ok,err]=ordinary_batch(R(rows,:),m,tolerance);
    Q(rows,:)=fallback;fast(rows)=ok;
    if nargout>2,residual(rows,:)=err;end
end

function [Q,fast,residual]=ordinary_batch(R,m,tolerance)
% Only comfortably interior, well-conditioned rows enter this array path.
% Original XYZ values remain available for the unchanged symbolic fallback.
    n=size(R,1);Q=nan(n,4);fast=false(n,1);residual=nan(n,2);
    scale=max(abs(R),[],2);valid=all(isfinite(R),2)&scale>0;
    rows=find(valid);if isempty(rows),return;end
    original=R(rows,:);scaled=original./scale(rows);r=scaled./sqrt(sum(scaled.^2,2));
    rho=hypot(r(:,2),r(:,3));beta=atan2(rho,r(:,1));
    outer=m.outer_deg*pi/180;inner=m.inner_deg*pi/180;go=outer-beta;gi=beta-inner;guard=32*eps;
    U2=sin((outer+beta)/2).*sin(go/2);V2=sin((beta+inner)/2).*sin(gi/2);
    lost=any(original~=0 & (r==0 | abs(r)<realmin),2);
    eligible=~lost & rho>0 & m.sin_product>0 & go>4*guard & gi>4*guard & ...
        isfinite(U2)&isfinite(V2)&U2>0&V2>0;
    ix=find(eligible);if isempty(ix),return;end
    U=sqrt(U2(ix));V=sqrt(V2(ix));budget=(guard*180/pi)*(1./U+1./V+1./rho(ix));
    keep=budget<=2.5e-9;ix=ix(keep);U=U(keep);V=V(keep);if isempty(ix),return;end
    eta=2*atan2d(U,V);delta=atan2d(2*U.*V,m.sin_alpha(1)*m.sin_difference+2*m.cos_alpha(1)*V2(ix));
    psi=atan2d(r(ix,3),r(ix,2));q=[psi+delta,-eta,psi-delta,eta];q(:,[1,3])=wrap(q(:,[1,3]));
    for j=[1,3]
        actual=q(:,j)+q(:,j+1);q(actual>180,j+1)=q(actual>180,j+1)-360;
        q(actual< -180,j+1)=q(actual< -180,j+1)+360;
    end
    sa=m.sin_alpha(1);sb=m.sin_alpha(2);ca=m.cos_alpha(1);cb=m.cos_alpha(2);rr=r(ix,:);
    % Evaluate the two independent forward chains in array columns. Preserve
    % the original arithmetic order and near-threshold row-wise fallback.
    u=q(:,[1,3]);v=q(:,[2,4]);su=sind(u);cu=cosd(u);sv=sind(v);cv=cosd(v);
    y=sa*cb+ca*sb.*cv;z=sb.*sv;
    fx=ca*cb-sa*sb.*cv;fy=cu.*y-su.*z;fz=su.*y+cu.*z;
    cx=fy.*rr(:,3)-fz.*rr(:,2);cy=fz.*rr(:,1)-fx.*rr(:,3);cz=fx.*rr(:,2)-fy.*rr(:,1);
    err=atan2d(sqrt((cx.^2+cy.^2)+cz.^2),(fx.*rr(:,1)+fy.*rr(:,2))+fz.*rr(:,3));
    % Borderline residual acceptance uses the original row-wise evaluator.
    keep=all(isfinite(q),2)&all(isfinite(err),2)&all(err<=min(tolerance/4,2.5e-11),2);
    dest=rows(ix(keep));Q(dest,:)=q(keep,:);fast(dest)=true;residual(dest,:)=err(keep,:);
end

function [q,status,dps,route]=solve(original,r,m)
    status='regular';dps=16;route='vector_half_angle';a=m.alpha_deg(1);b=m.alpha_deg(2);
    if original(2)==0 && original(3)==0
        if original(1)<0 || a~=b,unreachable();end
        q=[0,180,0,-180];status='center';route='vector_identity';return;
    end
    if original(1)==0 && m.outer_deg==90 && m.sum_roundoff_deg==0
        psi=atan2d(original(3),original(2));q=finish([psi,0,psi,0]);
        status='outer_boundary';route='vector_identity';return;
    end
    rho=hypot(r(2),r(3));beta=atan2(rho,r(1));outer=m.outer_deg*pi/180;inner=m.inner_deg*pi/180;
    go=outer-beta;gi=beta-inner;guard=32*eps;
    scaled_lost=any(original~=0 & (r==0 | abs(r)<realmin));
    uncertain=scaled_lost || rho==0 || m.sin_product==0 || abs(go)<=guard || abs(gi)<=guard;
    if ~uncertain && (go< -guard || gi< -guard),unreachable();end
    U2=sin((outer+beta)/2)*sin(go/2);V2=sin((beta+inner)/2)*sin(gi/2);
    if ~uncertain && U2>0 && V2>0 && all(isfinite([U2,V2]))
        U=sqrt(U2);V=sqrt(V2);
        if (guard*180/pi)*(1/U+1/V+1/rho)<=1e-8
            eta=2*atan2d(U,V);delta=atan2d(2*U*V,m.sin_alpha(1)*m.sin_difference+2*m.cos_alpha(1)*V2);
            psi=atan2d(r(3),r(2));q=finish([psi+delta,-eta,psi-delta,eta]);return;
        end
    end
    [handled,q,dps,native_status,native_route]=coneik_try_native_direction_precision(original,m);
    if handled,status=native_status;route=native_route;return;end
    [q,status,dps]=precise(original,m);route='vector_symbolic';
end

function [q,status,dps]=precise(r,m)
    if exist('sym','file')==0 || ~license('test','Symbolic_Toolbox')
        error('inverse_motion_analysis:PrecisionUnavailable','This difficult vector requires Symbolic Math Toolbox.');
    end
    prior_digits=digits;restore=onCleanup(@()digits(prior_digits)); %#ok<NASGU>
    exact=sym(r,'f');aa=sym(m.alpha_deg(1),'f');bb=sym(m.alpha_deg(2),'f');
    pi_exact=str2sym('pi');beta_exact=atan2(sqrt(exact(2)^2+exact(3)^2),exact(1));
    % Exact vector identities (including Pythagorean 45-degree targets) must
    % never be manufactured by a tolerance or rounded angle conversion.
    boundary='';
    if isAlways(beta_exact==(aa+bb)*pi_exact/180,'Unknown','false'),boundary='outer_boundary';
    elseif isAlways(beta_exact==abs(aa-bb)*pi_exact/180,'Unknown','false'),boundary='inner_boundary';end
    if ~isempty(boundary)
        digits(100);psi=double(vpa(atan2(exact(3),exact(2))*180/pi_exact,100));
        if strcmp(boundary,'outer_boundary'),q=finish([psi,0,psi,0]);
        else,if m.alpha_deg(1)<m.alpha_deg(2),psi=psi-180;end;q=finish([psi,-180,psi,180]);end
        status=boundary;dps=100;return;
    end
    previous='';previous_q=[];
    for dps=[100,200,400,800,1024]
        digits(dps);rad=vpa(pi_exact/180,dps);a=vpa(aa,dps)*rad;b=vpa(bb,dps)*rad;
        v=vpa(exact,dps);beta=atan2(sqrt(v(2)^2+v(3)^2),v(1));outer=a+b;inner=abs(a-b);
        go=outer-beta;gi=beta-inner;margin=vpa(sym(10)^(-dps+20),dps);
        if isAlways(abs(go)<=margin) || isAlways(abs(gi)<=margin),previous='';previous_q=[];continue;end
        if isAlways(go<0) || isAlways(gi<0)
            status='unreachable';raw=[];
        else
            U2=sin((outer+beta)/2)*sin(go/2);V2=sin((beta+inner)/2)*sin(gi/2);
            if ~isAlways(U2>0) || ~isAlways(V2>0),previous='';previous_q=[];continue;end
            U=sqrt(U2);V=sqrt(V2);eta=2*atan2(U,V)/rad;
            delta=atan2(2*U*V,sin(a)*sin(a-b)+2*cos(a)*V2)/rad;psi=atan2(v(3),v(2))/rad;
            raw=double([psi+delta,-eta,psi-delta,eta]);status='regular';
            if any(~isfinite(raw)),previous='';previous_q=[];continue;end
        end
        if strcmp(previous,status)
            if strcmp(status,'unreachable'),unreachable();end
            if max(abs(wrap(raw-previous_q)))<=1e-12,q=finish(raw);return;end
        end
        previous=status;previous_q=raw;
    end
    error('inverse_motion_analysis:NumericalUnresolved','Vector inverse could not be confirmed across precision levels.');
end

function q=finish(q)
    q([1,3])=wrap(q([1,3]));
    for j=[1,3]
        actual=q(j)+q(j+1);
        if actual>180,q(j+1)=q(j+1)-360;elseif actual< -180,q(j+1)=q(j+1)+360;end
    end
end
function v=wrap(v)
    v=rem(v,360);v(v>=180)=v(v>=180)-360;v(v< -180)=v(v< -180)+360;
end
function residual=assess(q,r,m)
    a=m.alpha_deg(1);b=m.alpha_deg(2);residual=zeros(1,2);
    for k=1:2
        u=q(2*k-1);v=q(2*k);y=sind(a)*cosd(b)+cosd(a)*sind(b)*cosd(v);z=sind(b)*sind(v);
        predicted=[cosd(a)*cosd(b)-sind(a)*sind(b)*cosd(v),cosd(u)*y-sind(u)*z,sind(u)*y+cosd(u)*z];
        residual(k)=atan2d(norm(cross(predicted,r)),dot(predicted,r));
    end
end
function unreachable()
    error('inverse_motion_analysis:Unreachable','Target direction is outside the model reachable annulus.');
end

function [q,residual,accepted]=ordinary_single(original,m,tolerance)
% Conservative success-only accelerator. Classification and difficult inputs
% remain the responsibility of the original engine, using original doubles.
    q=[];residual=[];accepted=false;scale=max(abs(original));
    if ~all(isfinite(original)) || scale==0,return;end
    scaled=original/scale;r=scaled/sqrt(sum(scaled.^2));
    rho=hypot(r(2),r(3));beta=atan2(rho,r(1));
    outer=m.outer_deg*pi/180;inner=m.inner_deg*pi/180;go=outer-beta;gi=beta-inner;guard=32*eps;
    if any(original~=0 & (r==0 | abs(r)<realmin)) || rho==0 || m.sin_product==0 || go<=4*guard || gi<=4*guard,return;end
    U2=sin((outer+beta)/2)*sin(go/2);V2=sin((beta+inner)/2)*sin(gi/2);
    if ~(U2>0 && V2>0 && isfinite(U2) && isfinite(V2)),return;end
    U=sqrt(U2);V=sqrt(V2);
    if (guard*180/pi)*(1/U+1/V+1/rho)>2.5e-9,return;end
    eta=2*atan2d(U,V);delta=atan2d(2*U*V,m.sin_alpha(1)*m.sin_difference+2*m.cos_alpha(1)*V2);
    psi=atan2d(r(3),r(2));q=finish([psi+delta,-eta,psi-delta,eta]);
    sa=m.sin_alpha(1);sb=m.sin_alpha(2);ca=m.cos_alpha(1);cb=m.cos_alpha(2);
    residual=zeros(1,2);
    for branch=1:2
        u=q(2*branch-1);v=q(2*branch);su=sind(u);cu=cosd(u);sv=sind(v);cv=cosd(v);
        y=sa*cb+ca*sb*cv;z=sb*sv;
        fx=ca*cb-sa*sb*cv;fy=cu*y-su*z;fz=su*y+cu*z;
        cx=fy*r(3)-fz*r(2);cy=fz*r(1)-fx*r(3);cz=fx*r(2)-fy*r(1);
        residual(branch)=atan2d(sqrt(cx*cx+cy*cy+cz*cz),fx*r(1)+fy*r(2)+fz*r(3));
    end
    % Near the acceptance threshold the original evaluator remains decisive.
    accepted=all(isfinite(q)) && all(isfinite(residual)) && all(residual<=min(tolerance/4,2.5e-11));
end

