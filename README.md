# ConeIK

ConeIK implements inverse kinematics for two-stage conical beam-steering mechanisms. It recovers candidate joint angles from beam directions and supports batch inversion of a trajectory.

Companion code for *ConeIK: Accurate joint-angle recovery and efficient batch inverse kinematics for two-stage conical beam-steering mechanisms*.

## Requirements

The MATLAB solver and synthetic example were tested with MATLAB R2024b. The Symbolic Math Toolbox is used by the variable-precision fallback for difficult inputs. The `code/native/` directory includes optional C++ MEX backends along with MATLAB build scripts in `code/`. Building them requires a MATLAB-supported C++ compiler. The two precision MEX targets also require Boost headers. Pass the Boost include root to the corresponding build function or set `BOOST_ROOT`.

## Quick start

Run the one-click entry from the repository root in MATLAB:

```matlab
main
```

`main` reads a 3-column trajectory txt (`time_s pitch_deg yaw_deg`, tabs/spaces/commas, one optional header line), solves every record, differentiates the recovered joint angles, and writes the result txt named in the `EDIT HERE` block at the top of `main.m` with the columns `time_s theta1_deg theta2_deg dtheta1_deg_s dtheta2_deg_s ddtheta1_deg_s2 ddtheta2_deg_s2`. Without arguments it runs the included 10,001-row synthetic trajectory (10 s at 1 kHz) in `data/`; edit `INPUT_FILE` in that block to point at your own txt. Custom calls: `main('my_trajectory.txt')`, `main(in, out, window)` also work. The automatic entry tries the GPU for batches of at least 4,096 rows and uses the CPU for smaller batches or when GPU execution is unavailable.

The underlying API, for reference:

```matlab
addpath('code');
model = coneik_model(20, 20);
Q = coneik_inverse_direction([1, 0.1, 0.1], model);

solve = coneik_prepare_direction(model);
Q_batch = solve([1, 0.1, 0.1; 1, -0.1, 0.1]);
```

Each row of the batch input contains a Cartesian beam direction. The four columns of the result contain two candidate joint-angle pairs in degrees. The original synthetic check remains available as `tests/demo_batch_inverse.m`.

## GPU batch

`code/coneik_prepare_direction_gpu.m` and `code/coneik_inverse_batch_gpu.m` explicitly request GPU execution for any batch size. They require the Parallel Computing Toolbox and a usable GPU device; devices newer than the bundled CUDA libraries are handled through driver forward compatibility, enabled automatically. The half-angle pipeline runs as one fused `arrayfun` kernel per batch. Rows the kernel declines (exact identities, boundaries, guarded or difficult inputs) fall back to the CPU row engine.

```matlab
solve_gpu = coneik_prepare_direction_gpu(model);
Q = solve_gpu(R);                          % N-by-3 directions
Q = coneik_inverse_batch_gpu(gamma, model); % N-by-2 [pitch, yaw] degrees
```

On a GeForce RTX 5050 against the vectorized CPU core (`tests/bench_gpu_batch.m`), the GPU becomes faster from a few thousand rows per batch and reaches 3.4x at 524,288 rows (0.105 vs 0.355 microseconds per target); candidate angles agree with the CPU core to 2e-13 degrees. Small batches stay on the CPU path.

[Recorded RTX 5050 synthetic benchmark](tests/results/bench_gpu_batch_rtx5050_public.txt).

## Data scope

This repository provides the solver and a synthetic demonstration. Access to and redistribution of the measured device-command data require permission from the data holder.
