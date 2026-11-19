function coneik_assert_gpu()
%CONEIK_ASSERT_GPU Fail fast when no usable GPU device is present.
% The GPU batch entry points call this once at prepare time so a missing
% device is reported before any target reaches the accelerator.
    try
        % Devices newer than the bundled CUDA libraries still work through
        % driver forward compatibility; enabling it twice is a no-op.
        parallel.gpu.enableCUDAForwardCompatibility(true);
    catch
    end
    try
        gpuDevice;
    catch me
        error('coneik:gpu:NoGPU', ...
            'A usable GPU device is required by the GPU batch entry points (%s).',me.message);
    end
end
