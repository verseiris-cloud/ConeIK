function coneik_position_selection_input(target,reference_q)
    if ~(isnumeric(target)&&isreal(target)&&isvector(target)&&~issparse(target)&&any(numel(target)==[2,3]))
        error('coneik_select_inverse_position:InvalidInput','target must be a full real 2- or 3-vector.');
    end
    if ~(isnumeric(reference_q)&&isreal(reference_q)&&isvector(reference_q)&&numel(reference_q)==2&&all(isfinite(reference_q(:))))
        error('coneik_select_inverse_position:InvalidReference','reference_q must be a finite real 2-vector.');
    end
end
