# ConeIK

ConeIK implements inverse kinematics for two-stage conical beam-steering mechanisms. It recovers candidate joint angles from beam directions and supports batch inversion of a trajectory.

Companion code for *ConeIK: Accurate joint-angle recovery and efficient batch inverse kinematics for two-stage conical beam-steering mechanisms*.

## Requirements

The MATLAB solver and synthetic example were tested with MATLAB R2024b. The Symbolic Math Toolbox is used by the variable-precision fallback for difficult inputs. The `code/native/` directory includes optional C++ MEX backends along with MATLAB build scripts in `code/`. Building them requires a MATLAB-supported C++ compiler. The two precision MEX targets also require Boost headers. Pass the Boost include root to the corresponding build function or set `BOOST_ROOT`.

## Quick start

From the repository root in MATLAB, run the synthetic batch example:

```matlab
addpath('tests');
demo_batch_inverse;
```

The example uses the 101-row synthetic trajectory in `data/sample_synthetic.dat`. Its three columns are time in seconds, pitch in degrees, and yaw in degrees. The results include the recovered joint angles and their velocities and accelerations.

To solve one beam direction and prepare a solver for a batch:

```matlab
addpath('code');
model = coneik_model(20, 20);
Q = coneik_inverse_direction([1, 0.1, 0.1], model);

solve = coneik_prepare_direction(model);
Q_batch = solve([1, 0.1, 0.1; 1, -0.1, 0.1]);
```

Each row of the batch input contains a Cartesian beam direction. The four columns of the result contain two candidate joint-angle pairs in degrees.

## Data scope

This repository provides the solver and a synthetic demonstration. Access to and redistribution of the measured device-command data require permission from the data holder.
