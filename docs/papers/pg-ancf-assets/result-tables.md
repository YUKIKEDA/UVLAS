| Group | finest n | DOFs | max energy | max position | max thin ratio* |
|---|---:|---:|---:|---:|---:|
| G1 | 4 | 1305 | 8.193649e-07 | 1.022500e-10 | 1.360525 |
| G2 | 4 | 1530 | 4.410056e-08 | 1.117758e-10 | 1.000000 |
| G3 | 4 | 1365 | 8.196357e-07 | 1.101195e-10 | 1.003037 |
| G4 | 4 | 1545 | 4.343268e-08 | 1.076631e-10 | 1.000000 |
| G5 | 2 | 708 | 8.035147e-06 | 5.466687e-08 | 0.829496 |
| G6 | 2 | 711 | 6.156903e-06 | 5.386094e-08 | 0.783333 |
| G7 | 2 | 714 | 9.563215e-06 | 5.412052e-08 | 0.853473 |

*Thin ratio is the maximum original per-RHS energy/position ratio over BOTH saved resolutions, with the original 1e-12 floor. It is not an asymptotic convergence bound.

| Domain | n | h | energy | position | original stationarity | load agreement |
|---|---:|---:|---:|---:|---:|---:|
| rectangle | 1 | 1e-02 | 7.361892e-05 | 9.918454e-04 | 6.313e-16 | 3.168e-13 |
| rectangle | 1 | 1e-04 | 7.544520e-07 | 6.974498e-05 | 3.341e-14 | 3.865e-10 |
| rectangle | 2 | 1e-02 | 1.442369e-05 | 4.927876e-05 | 7.011e-16 | 6.411e-13 |
| rectangle | 2 | 1e-04 | 1.478420e-07 | 1.072600e-05 | 3.551e-14 | 1.640e-09 |
| taper | 1 | 1e-02 | 7.626405e-05 | 9.885157e-04 | 5.814e-16 | 3.528e-13 |
| taper | 1 | 1e-04 | 7.781351e-07 | 8.462335e-05 | 1.392e-14 | 6.490e-10 |
| taper | 2 | 1e-02 | 1.582504e-05 | 5.650556e-05 | 6.324e-16 | 1.245e-12 |
| taper | 2 | 1e-04 | 1.601476e-07 | 1.425346e-05 | 2.691e-14 | 2.095e-09 |
