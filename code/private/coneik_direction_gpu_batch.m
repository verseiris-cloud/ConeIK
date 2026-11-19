function [Q,fast,residual]=coneik_direction_gpu_batch(R,m,tolerance)
%CONEIK_DIRECTION_GPU_BATCH GPU twin of the ordinary_batch array core.
% Only comfortably interior, well-conditioned rows are accepted here. Every
% declined row stays NaN and remains the caller's responsibility; the CPU
% row engine keeps identities, exact boundaries and the precision ladder.
% Thresholds and acceptance tests replicate ordinary_batch in
% coneik_direction_resolved: the same 32*eps guard, the same 2.5e-9 joint
% budget, and the same min(tolerance/4,2.5e-11) residual gate.
% The whole per-row pipeline runs as ONE fused arrayfun kernel, so a batch
% costs a single kernel launch plus the input/output transfers instead of
% one launch per elementwise operation. Degree trigonometry is expressed
% through radians; eps and realmin appear as literal constants because the
% kernel form of arrayfun does not call functions.
    n=size(R,1);Q=nan(n,4);fast=false(n,1);residual=nan(n,2);
    if n==0,return;end
    gate=min(tolerance/4,2.5e-11);
    [q1,q2,q3,q4,e1,e2]=arrayfun(@coneik_gpu_kernel, ...
        gpuArray(R(:,1)),gpuArray(R(:,2)),gpuArray(R(:,3)), ...
        m.outer_deg*pi/180,m.inner_deg*pi/180, ...
        m.sin_alpha(1),m.sin_alpha(2),m.cos_alpha(1),m.cos_alpha(2), ...
        m.sin_product,m.sin_difference,gate);
    q1=gather(q1);q2=gather(q2);q3=gather(q3);q4=gather(q4);
    e1=gather(e1);e2=gather(e2);
    fast=isfinite(q1)&isfinite(q2)&isfinite(q3)&isfinite(q4);
    Q(fast,:)=[q1(fast),q2(fast),q3(fast),q4(fast)];
    residual(fast,:)=[e1(fast),e2(fast)];
end

function [q1,q2,q3,q4,e1,e2]=coneik_gpu_kernel(x1,x2,x3,outer,inner,sa,sb,ca,cb,sp,sd,gate)
% Scalar per-row kernel compiled to CUDA by gpuArray arrayfun.
    guard=32*2.220446049250313e-16;
    realmin_=2.225073858507201e-308;
    deg=180/pi;
    q1=NaN;q2=NaN;q3=NaN;q4=NaN;e1=NaN;e2=NaN;
    scale=max(abs(x1),max(abs(x2),abs(x3)));
    if scale>0
        r1=x1/scale;r2=x2/scale;r3=x3/scale;
        nn=r1*r1+r2*r2+r3*r3;
        r1=r1/sqrt(nn);r2=r2/sqrt(nn);r3=r3/sqrt(nn);
        rho=hypot(r2,r3);beta=atan2(rho,r1);
        go=outer-beta;gi=beta-inner;
        U2=sin((outer+beta)/2)*sin(go/2);V2=sin((beta+inner)/2)*sin(gi/2);
        lost=(x1~=0&&(r1==0||abs(r1)<realmin_))||(x2~=0&&(r2==0||abs(r2)<realmin_))||(x3~=0&&(r3==0||abs(r3)<realmin_));
        if ~lost&&rho>0&&sp>0&&go>4*guard&&gi>4*guard&&U2>0&&V2>0
            U=sqrt(U2);V=sqrt(V2);
            if (guard*deg)*(1/U+1/V+1/rho)<=2.5e-9
                eta=2*atan2(U,V)*deg;
                delta=atan2(2*U*V,sa*sd+2*ca*V2)*deg;
                psi=atan2(r3,r2)*deg;
                a1=psi+delta;b1=-eta;a3=psi-delta;b3=eta;
                a1=rem(a1,360);a1=a1-360*(a1>=180)+360*(a1<-180);
                a3=rem(a3,360);a3=a3-360*(a3>=180)+360*(a3<-180);
                t1=a1+b1;if t1>180,b1=b1-360;elseif t1<-180,b1=b1+360;end
                t3=a3+b3;if t3>180,b3=b3-360;elseif t3<-180,b3=b3+360;end
                % Forward-chain residual for both branches, array form.
                rad=pi/180;
                y1=sa*cb+ca*sb*cos(b1*rad);z1=sb*sin(b1*rad);
                f1=ca*cb-sa*sb*cos(b1*rad);
                g1=cos(a1*rad)*y1-sin(a1*rad)*z1;h1=sin(a1*rad)*y1+cos(a1*rad)*z1;
                c1=g1*r3-h1*r2;c2=h1*r1-f1*r3;c3=f1*r2-g1*r1;
                ee1=atan2(sqrt(c1*c1+c2*c2+c3*c3),f1*r1+g1*r2+h1*r3)*deg;
                y3=sa*cb+ca*sb*cos(b3*rad);z3=sb*sin(b3*rad);
                f3=ca*cb-sa*sb*cos(b3*rad);
                g3=cos(a3*rad)*y3-sin(a3*rad)*z3;h3=sin(a3*rad)*y3+cos(a3*rad)*z3;
                d1=g3*r3-h3*r2;d2=h3*r1-f3*r3;d3=f3*r2-g3*r1;
                ee3=atan2(sqrt(d1*d1+d2*d2+d3*d3),f3*r1+g3*r2+h3*r3)*deg;
                if a1==a1&&b1==b1&&a3==a3&&b3==b3&&ee1==ee1&&ee3==ee3&&ee1<=gate&&ee3<=gate
                    q1=a1;q2=b1;q3=a3;q4=b3;e1=ee1;e2=ee3;
                end
            end
        end
    end
end
