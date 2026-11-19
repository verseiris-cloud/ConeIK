function [handled,q,dps]=coneik_try_native_precision(g1,g2,alpha)
% Optional precision accelerator; unresolved rows retain the original ladder.
    handled=false;q=[];dps=NaN;
    precision_backend=coneik_native_backend('precision');
    if isempty(precision_backend),return;end
    [row,code,precision]=precision_backend([g1,g2],alpha);
    if code==0
        q=reshape(row,2,2)';dps=precision;handled=true;
    elseif code==-1
        error('inverse_motion_analysis:Unreachable',coneik_precision_unreachable_message(alpha));
    elseif code~=-2
        error('inverse_motion_analysis:NumericalUnresolved','Unexpected native precision status.');
    end
end
