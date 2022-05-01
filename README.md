# Sparse stochastic matrix factorization


Code accompanying the paper "A Row-Wise Update Algorithm for Sparse Stochastic Matrix Factorization", Arxiv version available at:
https://arxiv.org/abs/2110.10412

README:
======================
Files marked 'xxx.mat' contain the numerical experiments in the paper.
The remaining files are standalone with no installation required.

1. Main algorithms

comparemain.m--Compare the effects of PALMSSMF, ALG24SSMF, and ARTMsparse algorithms.
PALMSSMF.m--Alg. 2.2 in the paper.
ALG24SSMF.m--Alg. 2.4 in the paper.

2. Auxiliary algorithms
SimplexProj.m--The algorithm of projecting the vector to the simplex corresponds to the Algorithm 3.3 in the paper.

3.data

V.mat--Corresponding to the statistical data of articles related to COVID-19, corresponding to example 4.5 in the article.
W0.mat, H0.mat--nitial point for example 4.5.
AW.mat,AH.mat--The result data obtained from the operation of algorithm ARTMsparse for example 4.5.
Atime.mat, Aerror.mat--Time and Hellinger distance data of the process of algorithm ARTMsparse for example 4.5.
TW.mat,TH.mat--The result data obtained from the operation of PALMSSMF for example 4.5.
Ttime.mat, Terror.mat--Time and Hellinger distance data of the process of PALMSSMF for example 4.5.
PW.mat,PH.mat--The result data obtained from the operation of ALG24SSMF for example 4.5.
Ptime.mat, Perror.mat--Time and Hellinger distance data of the process of ALG24SSMF for example 4.5.

4.figures
coviddata.eps coviddata.fig--figure of the comparison results of PALMSSMF, ALG24SSMF, and ARTMsparse algorithms
