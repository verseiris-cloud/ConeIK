function [handled,q]=coneik_unequal_axis(g1,g2,model)
% Restricted unequal-cone axial path, matching native/coneik_core.hpp.
% Exact degree gaps, no boundary snapping; unresolved/outside go to fallback.
handled=false;q=[];a=model.alpha_deg(1);b=model.alpha_deg(2);
if a==b || min(a,b)<1 || max(a,b)>80 || abs(a-b)<.25 || ...
        model.sum_roundoff_deg~=0 || model.difference_roundoff_deg~=0,return;end
if g1==0,beta=abs(g2);psi=-180*(g2<0);
elseif g2==0,beta=abs(g1);psi=-90*sign(g1);
elseif g1==-180 && abs(g2)>=90,beta=180-abs(g2);psi=-180*(g2>0);
elseif g2==-180 && abs(g1)>=90,beta=180-abs(g1);psi=-90*sign(g1);
else,return;end
outer=model.outer_deg;inner=model.inner_deg;
if ~((beta>=.9*outer && beta<=outer)||(beta>=inner && beta<=1.1*inner)),return;end
go=outer-beta;gi=beta-inner;if ~(go>0 && gi>0),return;end
U2=sind((outer+beta)/2)*sind(go/2);V2=sind((beta+inner)/2)*sind(gi/2);
if ~(U2>0&&V2>0),return;end
U=sqrt(U2);V=sqrt(V2);t1=sind(a)*sind(a-b);t2=2*cosd(a)*V2;
A=t1+t2;B=2*U*V;radius=hypot(A,B);
budget=128*eps*(180/pi)*(1+(abs(t1)+abs(t2)+B)/radius);
if ~(radius>0&&budget<1e-10),return;end
eta=2*atan2(U,V)*(180/pi);delta=atan2(B,A)*(180/pi);
q=[psi+delta,-eta;psi-delta,eta];
if any(abs(abs(sum(q,2))-180)<=16384*eps),q=[];return;end
handled=all(isfinite(q(:)));
end
