function [q11,q21,q12,q22,info] = coneik_inverse_position_resolved(gamma1,gamma2,model,allow_native_precision)
% Private arithmetic; MODEL was fully resolved by a parent-directory entry.
    if nargin<4,allow_native_precision=true;end
    if ~valid_scalar(gamma1) || ~valid_scalar(gamma2)
        error('inverse_motion_analysis:InvalidInput', ...
            'gamma1 and gamma2 must be finite real numeric scalars in degrees.');
    end
    g1 = wrap_degrees(double(gamma1)); g2 = wrap_degrees(double(gamma2));
    need_info=nargout>4;
    if need_info
        info = struct('status','regular','route','fast_half_general', ...
            'used_fallback',false,'precision_dps',16,'branch_count',2, ...
            'model_key',model.key);
    end
    [identified,q,status] = exact_identity(g1,g2,model);
    if identified
        if need_info
            info.status = status; info.route = 'identity';
            if strcmp(status,'center'), info.branch_count = Inf;
            else, info.branch_count = 1; end
        end
    else
        p = g1*(pi/180); y = g2*(pi/180);
        if abs(g1) <= 90 && abs(g2) <= 90
            cp = cos(p); x = cp*cos(y); ry = cp*sin(y); rz = -sin(p);
        else
            cp = cosd(g1); x = cp*cosd(g2); ry = cp*sind(g2); rz = -sind(g1);
        end
        rho = hypot(ry,rz); beta = atan2(rho,x);
        outer = model.outer_deg*(pi/180);
        inner = model.inner_deg*(pi/180);
        go = outer-beta; gi = beta-inner;
        % Difference-of-squares in factored half-angle form. Trigonometric
        % products retain tiny beta when direct 1-x would round to zero.
        U2 = sin((outer+beta)/2)*sin((outer-beta)/2);
        V2 = sin((beta+inner)/2)*sin((beta-inner)/2);
        guard = 32*eps; joint_budget_deg = 1e-8; reason = '';
        if any(abs([g1,g2,model.alpha_deg]) > 0 & abs([g1,g2,model.alpha_deg]) < realmin*180/pi)
            reason = 'subnormal_input';
        elseif rho == 0 || model.sin_product == 0
            reason = 'underflow';
        elseif abs(go) <= guard || abs(gi) <= guard
            reason = 'gap_guard';
        elseif go < -guard || gi < -guard
            error('inverse_motion_analysis:Unreachable', ...
                'Target direction lies outside the model annulus |alpha1-alpha2| <= beta <= alpha1+alpha2.');
        elseif ~all(isfinite([U2,V2])) || U2 <= 0 || V2 <= 0
            reason = 'radicand';
        elseif (180/pi)*guard*(1/sqrt(U2)+1/sqrt(V2)+1/rho) > joint_budget_deg
            reason = 'joint_budget';
        end
        if isempty(reason)
            U = sqrt(U2); V = sqrt(V2);
            eta = (180/pi)*2*atan2(U,V);
            A = model.sin_alpha(1)*model.sin_difference + 2*model.cos_alpha(1)*V2;
            B = 2*U*V;
            delta = (180/pi)*atan2(B,A); psi = (180/pi)*atan2(rz,ry);
            q = [psi+delta,-eta;psi-delta,eta];
        else
            handled=false;
            if model.alpha_deg(1)==model.alpha_deg(2)
                [handled,q]=coneik_ikgeo_axis(g1,g2,model.alpha_deg(1),false);
            else
                [handled,q]=coneik_unequal_axis(g1,g2,model);
            end
            if handled
                if need_info
                    info.route='ikgeo_sp4_sp1_axis';
                    if model.alpha_deg(1)~=model.alpha_deg(2),info.route='factored_unequal_axis';end
                end
            else
                [q,dps,native] = precise_inverse(g1,g2,model,allow_native_precision);
                if need_info
                    info.route = reason;if native,info.route='cpp_precision';end
                    info.used_fallback = true; info.precision_dps = dps;
                end
            end
        end
    end
    if ~all(isfinite(q(:)))
        error('inverse_motion_analysis:NumericalUnresolved','Nonfinite inverse output.');
    end
    q(:,1) = wrap_degrees(q(:,1));
    for k = 1:2
        if q(k,1)+q(k,2) > 180, q(k,2) = q(k,2)-360;
        elseif q(k,1)+q(k,2) < -180, q(k,2) = q(k,2)+360; end
    end
    q11=q(1,1); q21=q(1,2); q12=q(2,1); q22=q(2,2);
end

function [identified,q,status] = exact_identity(g1,g2,model)
    identified = false; q = []; status = '';
    a = model.alpha_deg(1); b = model.alpha_deg(2);
    center = (g1 == 0 && g2 == 0) || (g1 == -180 && g2 == -180);
    if center
        if a ~= b
            error('inverse_motion_analysis:Unreachable','The center direction is unreachable for unequal cone angles.');
        end
        identified = true; q = [0,180;0,-180]; status = 'center'; return;
    end
    % Represent beta as an exact sum of two binary64 degree values. Compare
    % expansions, not beta==fl(a+b): the latter can invent a boundary.
    if abs(g1) == 90
        bh=90; bl=0; psi=-90*sign(g1);
    elseif abs(g2) == 90
        bh=90; bl=0; psi=atan2d(-sind(g1),cosd(g1)*sind(g2));
    elseif g1 == 0
        bh=abs(g2); bl=0; psi=-180*(g2<0);
    elseif g2 == 0
        bh=abs(g1); bl=0; psi=-90*sign(g1);
    elseif g1 == -180
        [bh,bl]=two_sum(180,-abs(g2)); psi=-180*(g2>0);
    elseif g2 == -180
        [bh,bl]=two_sum(180,-abs(g1)); psi=-90*sign(g1);
    else
        return;
    end
    [oh,ol]=two_sum(a,b); [ih,il]=two_sum(max(a,b),-min(a,b));
    if bh==oh && bl==ol
        identified=true; q=[psi,0;psi,0]; status='outer_boundary';
    elseif bh==ih && bl==il
        if a < b, psi=psi-180; end
        identified=true; q=[psi,-180;psi,180]; status='inner_boundary';
    end
end

function [q,dps,native] = precise_inverse(g1,g2,model,allow_native_precision)
    native=false;
    if allow_native_precision
        [native,q,dps]=coneik_try_native_precision(g1,g2,model.alpha_deg);
        if native,return;end
    end
    if exist('sym','file') == 0 || ~license('test','Symbolic_Toolbox')
        error('inverse_motion_analysis:PrecisionUnavailable', ...
            'This difficult input requires Symbolic Math Toolbox for reliable classification and angles.');
    end
    saved_digits=digits; restore_digits=onCleanup(@() digits(saved_digits)); %#ok<NASGU>
    g1_exact=sym(g1,'f'); g2_exact=sym(g2,'f');
    a_exact=sym(model.alpha_deg(1),'f'); b_exact=sym(model.alpha_deg(2),'f');
    previous_state=''; previous_q=[];
    for dps=[100,200,400,800,1024]
        digits(dps); rad=vpa(str2sym('pi')/180,dps);
        a=vpa(a_exact,dps)*rad; b=vpa(b_exact,dps)*rad;
        p=vpa(g1_exact,dps)*rad; y=vpa(g2_exact,dps)*rad;
        x=cos(p)*cos(y); ry=cos(p)*sin(y); rz=-sin(p);
        beta=atan2(sqrt(ry*ry+rz*rz),x); inner=abs(a-b); outer=a+b;
        go=outer-beta; gi=beta-inner;
        margin=vpa(sym(10)^(-dps+20),dps);
        if isAlways(abs(go)<=margin) || isAlways(abs(gi)<=margin)
            previous_state=''; previous_q=[]; continue;
        end
        if isAlways(go<0) || isAlways(gi<0)
            state='unreachable'; q=[];
        else
            state='regular';
            U=sqrt(sin((outer+beta)/2)*sin((outer-beta)/2));
            V2=sin((beta+inner)/2)*sin((beta-inner)/2); V=sqrt(V2);
            eta=2*atan2(U,V)/rad;
            A=sin(a)*sin(a-b)+2*cos(a)*V2; B=2*U*V;
            delta=atan2(B,A)/rad; psi=atan2(rz,ry)/rad;
            q=double([psi+delta,-eta;psi-delta,eta]);
        end
        if strcmp(state,previous_state)
            if strcmp(state,'unreachable')
                error('inverse_motion_analysis:Unreachable', ...
                    'Target is outside the model annulus (checked at increasing precision).');
            end
            difference=abs(wrap_degrees(q-previous_q));
            if all(isfinite(q(:))) && max(difference(:))<=1e-12, return; end
        end
        previous_state=state; previous_q=q;
    end
    error('inverse_motion_analysis:NumericalUnresolved', ...
        'Classification or angles did not converge by 1024 decimal digits.');
end

function yes = valid_scalar(value)
    yes=isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value);
end

function [hi,lo] = two_sum(a,b)
    hi=a+b; z=hi-a; lo=(a-(hi-z))+(b-z);
end

function q = wrap_degrees(q)
    q=rem(q,360); q(q>=180)=q(q>=180)-360; q(q< -180)=q(q< -180)+360; q(q==0)=0;
end
