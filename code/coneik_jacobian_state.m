function state = coneik_jacobian_state(J)
%CONEIK_JACOBIAN_STATE Numerical rank and reusable SVD for a ConeIK Jacobian.
%   J is a finite real numeric 2x2 matrix, converted to double. Numerical
%   rank uses a machine-relative threshold; the 1e-8 singular-value-ratio
%   warning is separate and never changes rank. No damping is applied.
%   For a deficient J, V(:,active)*(U(:,active)'*target ./ s(active)) is the
%   minimum-norm retained-subspace solution (zero when active is empty).

    if ~(isnumeric(J) && isreal(J) && isequal(size(J),[2,2]) && ...
            all(isfinite(J(:))))
        error('coneik_jacobian_state:InvalidJacobian', ...
            'J must be a finite real numeric 2x2 matrix.');
    end
    J = double(full(J));
    [U,S,V] = svd(J);
    singular_values = diag(S);
    if ~all(isfinite(singular_values))
        error('coneik_jacobian_state:NumericalUnresolved', ...
            'Singular-value decomposition overflowed for the supplied matrix.');
    end
    sigma_max = singular_values(1);
    sigma_min = singular_values(2);
    if sigma_max == 0
        rank_tolerance = 0;
        active = false(2,1);
        condition_number = Inf;
        ratio = 0;
    else
        rank_tolerance = max(size(J))*eps(sigma_max);
        active = singular_values > rank_tolerance;
        ratio = sigma_min/sigma_max;
        if all(active), condition_number = sigma_max/sigma_min;
        else, condition_number = Inf; end
    end
    numerical_rank = sum(active);
    warning_ratio = 1e-8;
    state = struct('U',U,'s',singular_values,'V',V,'active',active,...
        'rank',numerical_rank,'condition',condition_number,...
        'near_singular',numerical_rank<2 || ratio<warning_ratio,...
        'rank_tolerance',rank_tolerance,'sigma_min',sigma_min,...
        'sigma_max',sigma_max,'near_singular_ratio_threshold',warning_ratio);
end
