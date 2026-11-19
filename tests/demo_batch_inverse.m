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

    assert(stats.n_total == 101, 'Unexpected synthetic sample length.');
    assert(stats.n_fail == 0, 'Synthetic sample contains unreachable rows.');
    assert(all(isfinite(result(:))), 'Batch result contains non-finite values.');

    max_step = max(max(abs(diff(result(:, 2:3), 1, 1))));
    fprintf('Synthetic batch check passed: %d/%d rows solved.\n', ...
        stats.n_ok, stats.n_total);
    fprintf('Maximum adjacent theta step: %.6f deg.\n', max_step);
    fprintf('First five output rows [t theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2]:\n');
    disp(result(1:min(5, size(result, 1)), :));
end
