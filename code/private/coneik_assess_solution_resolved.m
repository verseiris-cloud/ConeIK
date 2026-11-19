function assessment=coneik_assess_solution_resolved(gamma,Q,model,tolerance)
% Validated internal inputs. Evaluate both branches in array columns without
% N-by-3 temporary vectors. Keep the independent forward chain, scaled Gram
% calculation and factored small singular value used by the public assessor.
    valid_target=all(isfinite(gamma),2);g=rem(gamma,360);
    cp=cosd(g(:,1));tx=cp.*cosd(g(:,2));ty=cp.*sind(g(:,2));tz=-sind(g(:,1));
    beta=atan2d(hypot(ty,tz),tx);
    margin=min(beta-model.inner_deg,model.outer_deg-beta);margin(~valid_target)=NaN;
    sa=model.sin_alpha(1);ca=model.cos_alpha(1);sb=model.sin_alpha(2);cb=model.cos_alpha(2);
    u=rem(Q(:,1:2:end),360);v=rem(Q(:,2:2:end),360);
    su=sind(u);cu=cosd(u);sv=sind(v);cv=cosd(v);
    x=ca*cb-sa*sb.*cv;y0=sa*cb+ca*sb.*cv;z0=sb.*sv;
    y=cu.*y0-su.*z0;z=su.*y0+cu.*z0;
    rnorm=hypot(hypot(x,y),z);norm_error=abs(rnorm-1);
    available=valid_target & isfinite(u)&isfinite(v)&isfinite(rnorm)&rnorm>0&all(model.sin_alpha>0);
    cx=y.*tz-z.*ty;cy=z.*tx-x.*tz;cz=x.*ty-y.*tx;
    e=atan2d(hypot(hypot(cx,cy),cz),(x.*tx+y.*ty)+z.*tz);e(~available)=NaN;
    yp=-ca*sb.*sv;zp=sb.*cv;
    jx=sa*sb.*sv;jy=cu.*yp-su.*zp;jz=su.*yp+cu.*zp;
    scale=max(hypot(hypot(zeros(size(z)),z),y),hypot(hypot(jx,jy),jz));
    divisor=scale;divisor(divisor==0)=1;
    j1y=-z./divisor;j1z=y./divisor;
    j2x=jx./divisor;j2y=jy./divisor;j2z=jz./divisor;
    aa=j1y.^2+j1z.^2;bb=(j2x.^2+j2y.^2)+j2z.^2;ab=j1y.*j2y+j1z.*j2z;
    hi=scale.*sqrt(max(0,(aa+bb+hypot(aa-bb,2*ab))/2));
    lo=sa.*(sb./hi).*abs(sv);rank_loss=hi==0 | lo<=32*eps.*hi;lo(hi==0)=0;
    condition=hi./lo;gain=1./lo;condition(rank_loss)=Inf;gain(rank_loss)=Inf;
    hi(~available)=NaN;lo(~available)=NaN;condition(~available)=NaN;gain(~available)=NaN;norm_error(~available)=NaN;
    assessment=struct('direction_tolerance_deg',double(tolerance),'direction_error_deg',e,...
        'direction_pass',available & e<=tolerance,'unit_norm_error',norm_error,...
        'sigma_min',lo,'sigma_max',hi,'condition_number',condition,'inverse_gain',gain,...
        'numerically_singular',rank_loss & available,'evaluation_available',available,'valid_target',valid_target,...
        'boundary_margin_deg',margin,'boundary_margin_certified',false,'local_linearization_only',true,...
        'includes_model_or_measurement_uncertainty',false,'model_key',model.key);
end
