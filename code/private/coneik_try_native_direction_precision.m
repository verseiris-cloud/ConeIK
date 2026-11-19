function [handled,q,dps,status,route]=coneik_try_native_direction_precision(original,model)
% Only absence/unconfirmed classification delegates to the symbolic ladder.
handled=false;q=[];dps=NaN;status='regular';route='vector_cpp_precision';
fn=coneik_native_backend('direction_precision');
if isempty(fn),return;end
[q,code,dps]=fn(original,model.alpha_deg);
if code==0
    handled=true;q=reshape(q,1,4);
elseif code==1 || code==2
    handled=true;q=reshape(q,1,4);route='vector_cpp_exact_boundary';
    if code==1,status='outer_boundary';else,status='inner_boundary';end
elseif code==-1
    error('inverse_motion_analysis:Unreachable','Target direction is outside the model reachable annulus.');
elseif code~=-2
    error('inverse_motion_analysis:NumericalUnresolved','Native direction precision returned an unexpected classification.');
end
end
