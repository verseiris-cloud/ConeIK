function [score,detail] = coneik_state_selection_score(context,q,v,a)
%CONEIK_STATE_SELECTION_SCORE Internal soft ranking of already verified states.
% Context is validated by coneik_select_inverse_state; this is not a feasibility
% test. Physical q differences are not wrapped. Empty q omits position cost;
% fixed_q, when present, still contributes the constant position margin.
    parts=zeros(1,4);active=false(1,4);saturated=false;
    if ~isempty(q)&&~isempty(context.previous_q)
        [parts(1),sat]=component(q,context.previous_q,context.position_scale);
        saturated=saturated||sat;active(1)=true;
    end
    if ~isempty(context.preferred_velocity)
        [parts(2),sat]=component(v,context.preferred_velocity,context.velocity_scale);
        saturated=saturated||sat;active(2)=true;
    end
    if ~isempty(context.preferred_acceleration)
        [parts(3),sat]=component(a,context.preferred_acceleration,context.acceleration_scale);
        saturated=saturated||sat;active(3)=true;
    end
    qm=q;if isempty(qm)&&isfield(context,'fixed_q'),qm=context.fixed_q;end
    margins=ones(2,3);constrained=false(2,3);
    for k=1:2
        if ~isempty(qm)
            lo=context.bounds(k,1);hi=context.bounds(k,2);
            if isfinite(lo)&&isfinite(hi)
                constrained(k,1)=true;
                half=hi/2-lo/2;
                if half==0,margins(k,1)=0;
                else
                    margins(k,1)=min(delta(qm(k),lo,half),delta(hi,qm(k),half));
                end
            elseif isfinite(lo)
                constrained(k,1)=true;margins(k,1)=delta(qm(k),lo,360);
            elseif isfinite(hi)
                constrained(k,1)=true;margins(k,1)=delta(hi,qm(k),360);
            end
        end
        limits=[context.velocity_limits(k),context.acceleration_limits(k)];
        values=[v(k),a(k)];
        for j=1:2
            if isfinite(limits(j))
                constrained(k,j+1)=true;
                if limits(j)==0,margins(k,j+1)=0;
                else,margins(k,j+1)=1-abs(values(j))/limits(j);end
            end
        end
    end
    margins=max(0,min(1,margins));margin=min(margins(:));
    active(4)=any(constrained(:));parts(4)=1-margin;
    % Saturate only unrepresentable soft scores. Never change hard tests.
    cap=realmax/16;weighted=zeros(1,4);
    for j=1:4
        w=context.weights(j);
        if w==0||parts(j)==0,continue;end
        if parts(j)>cap/w,weighted(j)=cap;saturated=true;
        else,weighted(j)=parts(j)*w;end
    end
    score=sum(weighted);
    detail=struct('component_costs',parts,'weighted_costs',weighted,...
        'active_components',active,'normalized_margin',margin,...
        'axis_margins',margins,'constrained_axes',constrained,...
        'score',score,'score_saturated',saturated);
end

function [value,saturated]=component(x,reference,scale)
    d=delta(x,reference,scale);cap=sqrt(realmax/16);
    saturated=any(abs(d)>cap);d=max(-cap,min(cap,d));
    value=sum((d/sqrt(2)).^2);
end

function d=delta(x,y,scale)
    d=(x-y)./scale;
    overflow=~isfinite(x-y);
    if any(overflow),d(overflow)=x(overflow)./scale(overflow)-y(overflow)./scale(overflow);end
end
