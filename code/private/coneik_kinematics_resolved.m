function [r,gamma,J,b,detail] = coneik_kinematics_resolved(theta,dtheta,model)
% Private arithmetic kernel. Inputs and any explicit MODEL were validated
% by the current public adapter; MODEL=[] preserves implicit-default math.
    if ~isempty(model)
        if ~model.is_default
            switch nargout
                case {0,1}, r = general_kinematics(theta,dtheta,model);
                case 2, [r,gamma] = general_kinematics(theta,dtheta,model);
                case 3, [r,gamma,J] = general_kinematics(theta,dtheta,model);
                case 4, [r,gamma,J,b] = general_kinematics(theta,dtheta,model);
                otherwise, [r,gamma,J,b,detail] = general_kinematics(theta,dtheta,model);
            end
            return;
        end
    end

    % Reduce without adding 180 degrees: that addition loses tiny angles.
    theta = double(theta); dtheta = double(dtheta);
    q1 = rem(theta(1), 360); q2 = rem(theta(2), 360);
    s = sin(pi/9); c = cos(pi/9);
    su = sind(q1); cu = cosd(q1);
    sv = sind(q2); ch = cosd(q2/2);
    ch2 = ch*ch;
    % 1+cos(q2)=2*cos(q2/2)^2 avoids cancellation near the exact center.
    x = 1 - 2*s*s*ch2;
    y0 = 2*s*c*ch2; z0 = s*sv;
    y = cu*y0 - su*z0; z = su*y0 + cu*z0;
    r = [x; y; z];
    if nargout < 2, return; end

    cg = hypot(x, y);
    gamma = [-atan2d(z, cg); atan2d(y, x)];
    if nargout < 3, return; end

    cv = cosd(q2);
    xv = s*s*sv; yv0 = -s*c*sv; zv0 = s*cv;
    yv = cu*yv0 - su*zv0; zv = su*yv0 + cu*zv0;
    Jr = [0, xv; -z, yv; y, zv]; % d(direction)/d(theta_rad)
    D = x*x + y*y;
    J = [-Jr(3,:)/cg; (x*Jr(2,:) - y*Jr(1,:))/D];
    if nargout < 4, return; end

    v = dtheta(:);
    if all(v == 0) && nargout < 5
        b = [0; 0];
        return;
    end
    xvv = s*s*cv; yvv0 = -s*c*cv; zvv0 = -s*sv;
    yvv = cu*yvv0 - su*zvv0; zvv = su*yvv0 + cu*zvv0;
    % Three independent second-order coefficients, including the factor 2
    % in the mixed term. Derivatives below use radians; pi/180 converts the
    % final Hessian to degree output / degree-input squared.
    n1 = x*Jr(2,1) - y*Jr(1,1);
    n2 = x*Jr(2,2) - y*Jr(1,2);
    d1 = 2*(x*Jr(1,1) + y*Jr(2,1));
    d2 = 2*(x*Jr(1,2) + y*Jr(2,2));
    coefficients = [ ...
        z/cg-z*Jr(3,1)^2/cg^3, ...
        -2*yv/cg-2*z*Jr(3,1)*Jr(3,2)/cg^3, ...
        -zvv/cg-z*Jr(3,2)^2/cg^3; ...
        -x*y/D-n1*d1/D^2, ...
        -2*x*zv/D-(n1*d2+n2*d1)/D^2, ...
        (x*yvv-y*xvv)/D-n2*d2/D^2]*(pi/180);
    b = coefficients*[v(1)^2; v(1)*v(2); v(2)^2];
    if nargout > 4
        detail = struct('direction_jacobian_rad', Jr, ...
            'curvature_coefficients', coefficients, ...
            'curvature_coefficient_units', 'degree / degree^2', ...
            'angle_units', 'degrees', 'coordinate_system', 'nominal');
    end
end

function [r,gamma,J,b,detail] = general_kinematics(theta,dtheta,model)
    if any(model.sin_alpha==0)
        error('coneik_kinematics:NumericalUnresolved', ...
            'A positive cone angle underflows in binary64; direction and motion cannot be resolved by this numeric kernel.');
    end
    theta=double(theta); dtheta=double(dtheta);
    q1=rem(theta(1),360); q2=rem(theta(2),360);
    sa=model.sin_alpha(1); sb=model.sin_alpha(2); ca=model.cos_alpha(1);
    su=sind(q1); cu=cosd(q1); sv=sind(q2); cv=cosd(q2); ch=cosd(q2/2); ch2=ch*ch;
    % The inner-boundary representation retains the equal-angle center.
    x=model.cos_difference-2*sa*sb*ch2;
    y0=model.sin_difference+2*ca*sb*ch2; z0=sb*sv;
    if q2==0
        x=model.cos_sum; y0=model.sin_sum; z0=0;
    elseif abs(q2)==180
        x=model.cos_difference; y0=model.sin_difference; z0=0;
    end
    y=cu*y0-su*z0; z=su*y0+cu*z0; r=[x;y;z];
    if nargout<2, return; end
    cg=hypot(x,y);
    if cg==0
        error('coneik_kinematics:ChartSingularity', ...
            'Pitch/yaw is undefined at this direction pole; request direction only or use another output chart.');
    end
    gamma=[-atan2d(z,cg);atan2d(y,x)];
    if nargout<3, return; end
    xv=sa*sb*sv; yv0=-ca*sb*sv; zv0=sb*cv;
    yv=cu*yv0-su*zv0; zv=su*yv0+cu*zv0;
    Jr=[0,xv;-z,yv;y,zv]; D=x*x+y*y;
    J=[-Jr(3,:)/cg;(x*Jr(2,:)-y*Jr(1,:))/D];
    if nargout<4, return; end
    v=dtheta(:);
    if all(v==0)&&nargout<5, b=[0;0]; return; end
    xvv=sa*sb*cv; yvv0=-ca*sb*cv; zvv0=-sb*sv;
    yvv=cu*yvv0-su*zvv0; zvv=su*yvv0+cu*zvv0;
    n1=x*Jr(2,1)-y*Jr(1,1); n2=x*Jr(2,2)-y*Jr(1,2);
    d1=2*(x*Jr(1,1)+y*Jr(2,1)); d2=2*(x*Jr(1,2)+y*Jr(2,2));
    coefficients=[z/cg-z*Jr(3,1)^2/cg^3, ...
        -2*yv/cg-2*z*Jr(3,1)*Jr(3,2)/cg^3, ...
        -zvv/cg-z*Jr(3,2)^2/cg^3; ...
        -x*y/D-n1*d1/D^2, ...
        -2*x*zv/D-(n1*d2+n2*d1)/D^2, ...
        (x*yvv-y*xvv)/D-n2*d2/D^2]*(pi/180);
    b=coefficients*[v(1)^2;v(1)*v(2);v(2)^2];
    if nargout>4
        detail=struct('direction_jacobian_rad',Jr,'curvature_coefficients',coefficients, ...
            'curvature_coefficient_units','degree / degree^2','angle_units','degrees', ...
            'coordinate_system','nominal','model_key',model.key,'alpha_deg',model.alpha_deg);
    end
end
