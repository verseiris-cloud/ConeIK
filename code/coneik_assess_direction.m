function assessment=coneik_assess_direction(R,Q,model,options)
%CONEIK_ASSESS_DIRECTION Optional residual and local sensitivity checks for XYZ.
% A=CONEIK_ASSESS_DIRECTION(R,Q,MODEL,OPTIONS) assesses original N-by-3 beam
% directions against N-by-2 or N-by-4 NOMINAL joint angles in degrees. For
% one selected pair q (2-by-1), pass q.' as Q. A spatial point P must first
% be expressed as R=P-origin in the model coordinate frame.
%
% OPTIONS.direction_tolerance_deg: finite scalar in (0,180), default 1e-8.
% OPTIONS.direction_uncertainty_deg: optional scalar or N-by-1 assumed
% angular perturbation magnitude in [0,180]. Omitted means UNKNOWN, not 0.
% It is not an XYZ standard deviation or a calibrated measurement accuracy.
%
% Uses the original XYZ for the residual; no pitch/yaw conversion or inverse
% call. Positive rescaling of R leaves diagnostics unchanged up to rounding.
% Zero/nonfinite directions and nonfinite candidate pairs are unavailable.
% Normalization that loses or leaves subnormal nonzero components is flagged
% and not accepted as a diagnostic verdict. This double-precision assessor
% does not replace the inverse solver's multiprecision boundary decisions.
%
% Per-branch fields:
%   direction_error_deg, direction_pass: forward residual / tolerance test.
%   sigma_min, sigma_max: singular values of d(unit direction)/d(q_rad).
%   inverse_gain: 1/sigma_min; local joint 2-norm error / direction angle
%     error, in deg/deg. Inf denotes numerical rank loss, not zero error.
%   local_joint_error_estimate_deg: inverse_gain * supplied uncertainty.
%     A first-order worst-direction estimate, NOT a finite-error bound or
%     confidence interval. NaN if uncertainty is absent or residual fails.
%   local_estimate_available: regular, residual-passing, finite estimate.
% Per-target fields:
%   boundary_margin_deg: approximate angular distance to reachable cones.
%   uncertainty_reaches_boundary: supplied perturbation reaches that margin;
%     nonlinear effects/unreachable perturbations then need explicit study.
%   branch_separation_deg: periodic joint 2-norm distance for two supplied
%     pairs; NaN for one pair. This is not a proof of branch completeness.
%
% No branch is selected and no joint limits are checked here. The optional
% uncertainty describes direction only; geometry/calibration errors are not
% included. Small residual, low sensitivity and unique solution differ.
    if nargin<3,model=[];end
    if nargin<4,options=struct();end
    if ~(isnumeric(R)&&isreal(R)&&ismatrix(R)&&~issparse(R)&&size(R,2)==3)
        error('coneik_assess_direction:InvalidInput','R must be a full real N-by-3 numeric matrix.');
    end
    n=size(R,1);
    if ~(isnumeric(Q)&&isreal(Q)&&ismatrix(Q)&&~issparse(Q)&&size(Q,1)==n&&any(size(Q,2)==[2,4]))
        error('coneik_assess_direction:InvalidInput','Q must be a full real N-by-2 or N-by-4 matrix matching R rows.');
    end
    if ~(isstruct(options)&&isscalar(options))||...
            any(~ismember(fieldnames(options),{'direction_tolerance_deg','direction_uncertainty_deg'}))
        error('coneik_assess_direction:InvalidOptions','Unknown option or non-scalar options structure.');
    end
    tolerance=1e-8;
    if isfield(options,'direction_tolerance_deg'),tolerance=options.direction_tolerance_deg;end
    if ~(isnumeric(tolerance)&&isreal(tolerance)&&isscalar(tolerance)&&isfinite(tolerance)&&tolerance>0&&tolerance<180)
        error('coneik_assess_direction:InvalidTolerance','Tolerance must be a finite scalar in (0,180) degrees.');
    end
    uncertainty=nan(n,1);provided=isfield(options,'direction_uncertainty_deg');
    if provided
        uncertainty=options.direction_uncertainty_deg;
        if ~(isnumeric(uncertainty)&&isreal(uncertainty)&&~issparse(uncertainty)&&...
                (isscalar(uncertainty)||isequal(size(uncertainty),[n,1]))&&...
                all(isfinite(uncertainty(:)))&&all(uncertainty(:)>=0)&&all(uncertainty(:)<=180))
            error('coneik_assess_direction:InvalidUncertainty','Uncertainty must be a scalar or N-by-1 finite vector in [0,180] degrees.');
        end
        uncertainty=double(uncertainty);
        if isscalar(uncertainty),uncertainty=repmat(uncertainty,n,1);end
    end
    model=coneik_resolve_model(model);R=double(R);Q=double(Q);
    magnitude=max(abs(R),[],2);valid_target=all(isfinite(R),2)&magnitude>0;
    divisor=magnitude;divisor(~valid_target)=1;
    r=R./divisor;r=r./hypot(hypot(r(:,1),r(:,2)),r(:,3));
    limited=valid_target&any(R~=0 & (r==0 | abs(r)<realmin),2);
    target_available=valid_target&~limited;r(~target_available,:)=NaN;
    tx=r(:,1);ty=r(:,2);tz=r(:,3);
    beta=atan2d(hypot(ty,tz),tx);
    margin=min(beta-model.inner_deg,model.outer_deg-beta);margin(~target_available)=NaN;
    sa=model.sin_alpha(1);ca=model.cos_alpha(1);sb=model.sin_alpha(2);cb=model.cos_alpha(2);
    u=rem(Q(:,1:2:end),360);v=rem(Q(:,2:2:end),360);
    su=sind(u);cu=cosd(u);sv=sind(v);cv=cosd(v);
    x=ca*cb-sa*sb.*cv;y0=sa*cb+ca*sb.*cv;z0=sb.*sv;
    y=cu.*y0-su.*z0;z=su.*y0+cu.*z0;
    rnorm=hypot(hypot(x,y),z);norm_error=abs(rnorm-1);
    available=target_available&isfinite(u)&isfinite(v)&isfinite(rnorm)&rnorm>0&all(model.sin_alpha>0);
    cx=y.*tz-z.*ty;cy=z.*tx-x.*tz;cz=x.*ty-y.*tx;
    residual=atan2d(hypot(hypot(cx,cy),cz),(x.*tx+y.*ty)+z.*tz);
    residual(~available)=NaN;norm_error(~available)=NaN;
    yp=-ca*sb.*sv;zp=sb.*cv;
    jx=sa*sb.*sv;jy=cu.*yp-su.*zp;jz=su.*yp+cu.*zp;
    scale=max(hypot(z,y),hypot(hypot(jx,jy),jz));d=scale;d(d==0)=1;
    j1y=-z./d;j1z=y./d;j2x=jx./d;j2y=jy./d;j2z=jz./d;
    aa=j1y.^2+j1z.^2;bb=(j2x.^2+j2y.^2)+j2z.^2;ab=j1y.*j2y+j1z.*j2z;
    hi=scale.*sqrt(max(0,(aa+bb+hypot(aa-bb,2*ab))/2));
    % Product of singular values = sa*sb*abs(sin(q2)); avoid subtracting
    % almost equal eigenvalues of the Gram matrix near rank loss.
    lo=sa.*(sb./hi).*abs(sv);lo(hi==0)=0;
    sensitivity_available=available&isfinite(hi)&isfinite(lo);
    rank_loss=(hi==0|lo<=32*eps.*hi)&sensitivity_available;
    condition=hi./lo;gain=1./lo;condition(rank_loss)=Inf;gain(rank_loss)=Inf;
    hi(~sensitivity_available)=NaN;lo(~sensitivity_available)=NaN;
    condition(~sensitivity_available)=NaN;gain(~sensitivity_available)=NaN;
    pass=available&residual<=tolerance;
    estimate=gain.*uncertainty;
    if provided,estimate(rank_loss)=Inf;end
    estimate(~pass)=NaN;
    estimate_available=pass&sensitivity_available&~rank_loss&provided&isfinite(estimate);
    separation=nan(n,1);
    if size(Q,2)==4
        diff=rem(Q(:,1:2),360)-rem(Q(:,3:4),360);diff=mod(diff+180,360)-180;
        separation=hypot(diff(:,1),diff(:,2));separation(~all(isfinite(Q),2))=NaN;
    end
    assessment=struct('input_representation','direction_xyz','model_key',model.key,...
        'direction_tolerance_deg',double(tolerance),'direction_error_deg',residual,'direction_pass',pass,...
        'unit_norm_error',norm_error,'sigma_min',lo,'sigma_max',hi,'condition_number',condition,...
        'inverse_gain',gain,'numerically_singular',rank_loss,'evaluation_available',available,...
        'sensitivity_available',sensitivity_available,'valid_target',valid_target,...
        'normalization_precision_limited',limited,'boundary_margin_deg',margin,'boundary_margin_certified',false,...
        'branch_separation_deg',separation,'has_direction_uncertainty',provided,...
        'direction_uncertainty_deg',uncertainty,'local_joint_error_estimate_deg',estimate,...
        'local_estimate_available',estimate_available,...
        'uncertainty_reaches_boundary',target_available&provided&uncertainty>=margin,...
        'local_linearization_only',true,'includes_model_uncertainty',false);
end
