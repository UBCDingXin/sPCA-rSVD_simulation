# sPCA-rSVD_simulation
Codes to test the performance of sparse principal component analysis via regularized low rank matrix approximation (sPCA-rSVD) on synthetic datasets. (Shen, H., & Huang, J. Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation. Journal of multivariate analysis, 99(6), 1015-1034.)

"sPCA-rSVD.cpp" includes the necessary function to implement sPCA-rSVD and its cross-validation version.

"Simulation_Example1.R" tries to reproduce the results of Example 1 in the paper.

"Simulation_Robust.R" tests the robustness of sPCA-rSVD to random corruption.

"Simulation_SupervisedPCA.R" compares the performance of sPCA-rSVD, supervised PCA and Lasso on a sparse linear setting.
