How to use the Code...

**Corresponding paper:
Gautschi-type and implicit-explicit integrators for constrained wave equations 
by R. Altmann, B. Dörich, C. Zimmer

**Step 1: run experiments
--> run_kineticWave.m
this will run the Euler/CN/Gautschi schemes and compare it to a reference solution
code needs existing folders
- errors_PDAE_Euler
- errors_PDAE_CNIMEX
- errors_PDAE_Gautschi

options: 
- in line 12 you may change the spatial grid
- in line 14 you may change the temporal grid

**Step 2: plot results
--> run_convPlots.m
- plots the convergence results of one scheme (lines 12-14)
