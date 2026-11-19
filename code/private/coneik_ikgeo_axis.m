function [handled,q]=coneik_ikgeo_axis(g1,g2,a,preserve_legacy_error)
% SP4 followed by SP1, specialized to equal cones and an exact axial target.
% Only used near the outer boundary after the ordinary solver's guard fires.
% Original binary64 degrees retain the boundary gap; no acos, clamping,
% neighbourhood snapping, or least-squares solution is used.
    if nargin<4,preserve_legacy_error=true;end
    handled=false;q=[];direct=true;
    if a<1 || a>80,return;end
    if g1==0
        beta=abs(g2);psi=-180*(g2<0);
    elseif g2==0
        beta=abs(g1);psi=-90*sign(g1);
    elseif g1==-180 && abs(g2)>=90
        % Exact subtraction: both operands lie in [90,180] (Sterbenz).
        beta=180-abs(g2);psi=-180*(g2>0);direct=false;
    elseif g2==-180 && abs(g1)>=90
        beta=180-abs(g1);psi=-90*sign(g1);direct=false;
    else
        return;
    end
    outer=2*a;
    if beta<0.9*outer || beta>1.1*outer,return;end
    gap=outer-beta; % Exact subtraction by Sterbenz in this restricted range.
    if gap<0
        if preserve_legacy_error && direct
            error('inverse_motion_analysis:Unreachable', ...
                'Axial target is outside the equal-cone boundary (exact degree comparison).');
        end
        return; % New coverage only accepts positive gaps; keep original errors.
    end
    if gap==0,return;end % Exact identities belong to the existing identity route.
    % SP4 discriminant, factored before evaluation. U/sin(a)=sin(|q2|/2),
    % V/sin(a)=cos(q2/2); sqrt(norm_A^2-b^2) would lose the small gap.
    U=sqrt(sind((outer+beta)/2)*sind(gap/2));V=sind(beta/2);
    eta=2*atan2(U,V)*(180/pi);
    % SP1 about x: subtract the azimuth of the intermediate vector from psi.
    % Common factors cancel, leaving a stable atan2 even for tiny q2.
    delta=atan2(U,cosd(a)*V)*(180/pi);
    q=[psi+delta,-eta;psi-delta,eta];
    if ~(preserve_legacy_error && direct)
        % Preserve the old precision solver's periodic representative near
        % an actual-angle cut (notably negative yaw at a=60). Never snap it.
        if any(abs(abs(sum(q,2))-180)<=16384*eps),q=[];return;end
    end
    handled=all(isfinite(q(:)))&&U>0&&V>0;
end
