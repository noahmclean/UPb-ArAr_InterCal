%% Import data

addpath(genpath('../data/'))

datePairs = csvread('RenneEtAl2010_Table2.csv', 1, 0);
% From table supplied by Paul Jan 23rd
% Added 1 yr uncertainty for Vesuvius
% columns: age (Ma), age_sigma, R, r_sigma
% rows: same as in Renne et al 2010 table 2 
% (note row 1 / 79CE is calendar year, not U-Pb)

K40decay = Kdecay();
% data from Endt and Van der Leun compilation, recalculated by Min et al. (2000)
% [lambdaEC, 1s, lambdaBeta, 1s, rhoEC_Beta, lambdaTotal, 1s] - all unct. absolute

kappaFCs = FCsData();
% following Renne et al. 2010, K-Ar from FCs and through R_FCs-Std
% [kappaFCs, 1s, dkappaFCs_dr40K_K]


%% Replicate Kwon (2002) approach as described in Renne et al. 2010 eqns 

x0 = [K40decay(1,3) kappaFCs(1) datePairs(:,3)'];
