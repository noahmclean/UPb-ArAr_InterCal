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

kappaFCs = FCsData();