# Accelerated AltMin for SRPCP

## Overview

This project solves the Square-Root Principal Component Pursuit (SRPCP) problem

$$
\min_{L,S} \quad \lVert L\rVert_* + \lambda\lVert S\rVert_1 + \mu\lVert L+S-D\rVert_F,
$$

for $D\in\mathbb{R}^{m\times n}$, with $\lambda = 1/\sqrt{m}$ and $\mu = \sqrt{n/2}$, by
alternating minimization. The original implementation can be found at
[srPCP_code](https://github.com/MatOpt/srPCP_code).

This solver accelerates the sparse subproblem by replacing the full sort with an adaptive interval sort strategy. It outperforms the original solver in the tested large-scale, highly unbalanced cases. 

## Structure

```
src/
  altmin.m            outer AltMin loop
  updateS.m           sparse subproblem dispatcher
  updateS_full.m        full-sort solver (baseline)
  updateS_acc.m         interval-sort solver
  updateL.m           low-rank subproblem
  biquickselect.cpp   two-pivot selection, MEX source
  quickselect.cpp     selection, MEX source
tests/                experiment scripts
data/                 input video datasets (not included; see Data)
results/              generated figures and logs
```

## Build

- Developed with MATLAB R2023a.
- The accelerated solver uses two MEX files, built from `src/`. Note that the sources use C++17 parallel algorithms (Intel TBB on Linux).

  ```matlab
  cd src
  mex CXXFLAGS='$CXXFLAGS -std=c++17' biquickselect.cpp -ltbb
  mex CXXFLAGS='$CXXFLAGS -std=c++17' quickselect.cpp   -ltbb
  ```


## Usage

The following code shows a sample of running the algorithm developed in this project: 

```matlab
addpath(genpath("src"))

m = size(D, 1);
n = size(D, 2);
lambda = 1 / sqrt(m);
mu = sqrt(n / 2);

params = struct("Smod", "partialsort");
[L, S, res] = altmin(D, lambda, mu, params);
```

`params` fields:

| field     | meaning               | values / default                        |
|-----------|-----------------------|-----------------------------------------|
| `Lmod`    | low-rank solver mode  | `"partialsvd"` (default), `"fullsvd"`   |
| `Smod`    | sparse solver mode    | `"partialsort"` (default), `"fullsort"` |
| `tol`     | convergence tolerance | `1e-5`                                   |
| `maxiter` | maximum iterations    | `5000`                                   |
| `disp`    | per-iteration log     | `true`                                   |

`res` returns the status, final objective, total time, and per-subproblem time breakdown.

## Data

Video tests use the low-light raw videos from the DRV dataset (Chen et al., *Seeing Motion in the Dark*, ICCV 2019).

This repository does not contain the data due to its size. To obtain it, please run the following commands in the root path of the project:

```bash
mkdir -p data
wget -P ./data https://storage.googleapis.com/isl-datasets/DRV/VBM4D_rawRGB.zip
unzip ./data/VBM4D_rawRGB.zip -d ./data
rm ./data/VBM4D_rawRGB.zip
```
