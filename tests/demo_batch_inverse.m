function demo_batch_inverse()
%DEMO_BATCH_INVERSE Run batch inversion on the included synthetic sample.
% The sample is generated data and contains no measured H1 records.

    this_file = mfilename('fullpath');
    tests_dir = fileparts(this_file);
    package_root = fileparts(tests_dir);
    code_dir = fullfile(package_root, 'code');
    data_file = fullfile(package_root, 'data', 'sample_synthetic.dat');

    addpath(code_dir);

    [result, stats] = batch_inverse_dat(data_file, 7);

    assert(stats.n_total == 10001, 'Unexpected synthetic sample length.');
    assert(stats.n_fail == 0, 'Synthetic sample contains unreachable rows.');
    assert(all(isfinite(result(:))), 'Batch result contains non-finite values.');

    % Physical step size: unwrap first so a +/-360 degree representative
    % jump at the branch cut does not masquerade as a discontinuity.
    max_step = max(max(abs(rad2deg(diff(unwrap(deg2rad(result(:, 2:3))))))));
    fprintf('Synthetic batch check passed: %d/%d rows solved.\n', ...
        stats.n_ok, stats.n_total);
    fprintf('Maximum adjacent theta step: %.6f deg.\n', max_step);
    fprintf('First five output rows [t theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2]:\n');
    disp(result(1:min(5, size(result, 1)), :));
end
