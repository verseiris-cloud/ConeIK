function message=coneik_precision_unreachable_message(alpha)
% Preserve the symbolic branch's public error message.
    if all(alpha==[20,20])
        message='Target direction is outside beta = 40 degrees (checked at increasing precision).';
    else
        message='Target is outside the model annulus (checked at increasing precision).';
    end
end
