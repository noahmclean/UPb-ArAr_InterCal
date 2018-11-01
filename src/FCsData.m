function kappaFCs = FCsData()

%% Compilation of Fish Canyon Sanidine 40Ar*/40K (= kappaFCs) measurements
%  compile dnd convert reported data into mol/mol basis
%  keep track of 40K/K if possible.

%% Raw measurement info

% Steven et al. (1967)            
Steven = [11.3 NaN 4.68e-10 NaN 27.8 0.8/2 1.22e-4];
%        [K2O (%), 1s abs unct, Ar* (mol/g), 1s abs unct, ...
%                        date (Ma), 1s abs unct, 40K/K (g/g K)]

% Hurford and Hammerschmidt (1985)
HandH = [10.21 NaN (11.23+11.23)/2 NaN (27.93+28.04)/2 0.54/2 0.01167/100];
%       [%K, unct, 40Ar* (10^-6 cm^3 STP g^-1), date (Ma), unct, 40K/K (atom %)]
% Ar* is average of two measurements reported in Table 1. 
% Two ages are reported, one for each Ar* measurement, both ages reported with same unct.
% Unct in the (avg) age is taken as this reported uncertainty in each measurement

% Jourdan and Renne (2007) FCs vs. primary standards
JandR = zeros(4,6);
% HB3gr, K-Ar data from Turner et al. (1971)
JandR(1,:) = [1.2465 0.0030 7.09e-5 0.04e-5 0.019276 0.000028];
%            [wt% K, 1s abs unct, cm^3 STP/g, 1s abs, R_FCs, 1s abs unct]

% GA-1550, K ID data from Renne et al. (1998), Ar* from McDougall and Roksandic (1974)
JandR(2,:) = [7.626 0.016 NaN NaN 0.279803 0.00024];
% still looking for Ar* data

% NL-25 hornblende from amphibolite. 
% K from Schaefer and Schaefer (1977), Ar* from Alexander and Davis (1974)
JandR(3,:) = [0.3070 0.0008 7200e-8 40e-8 0.004730 0.000046];
%            [wt% K, 1s abs, cm^3 STP/g, 1s abs, R_FCs, 1s abs

% GHC-305 biotite from Sierran pluton
% K and Ar from Renne et al. (1998)
JandR(4,:) = [7.570 0.011 1.428e-9 0.004e-9 0.25947 0.00086];
%            [wt% K, 1s unct, mol Ar*/g, 1s unct, R_FCs, 1s abs


%% Calculating all data to same basis:

%% First, recreate Steven et al. (1967) FCs K-Ar date:

wtPctK2O_to_mol40KperGramOld=(1/100)*(1/(2*39.102+15.9994))*(2/1)*(Steven(7)*39.102/39.964);
% (g K20/100 g sample)*(mol K2O/grams K2O)*(mol K/mol K2O)*(mol 40K/mol K)
% note: converting g 40K/g K reported in Steven footnote to mol 40K/mol K using
% atomic weight of K from 1967 = 39.102, K2O atomic mass also using 1967 atomic weights.
% difference = 78 ppm bewteen old and current (2018) K2O and K atomic weights

Steven_mol40KperGramOld = Steven(1)*wtPctK2O_to_mol40KperGramOld;
Steven_kappaFCsOld = Steven(3)/Steven_mol40KperGramOld;

lambdaEC_Steven = 0.584e-10; %/yr, from footnote in Steven Table 2
lambdaBeta_Steven = 4.72e-10; %/yr, from footnote in Steven Table 2
lambdaTotal_Steven = lambdaEC_Steven + lambdaBeta_Steven;
StevensAgeOld = 1/lambdaTotal_Steven * ...
                           log(lambdaTotal_Steven/lambdaEC_Steven * Steven_kappaFCsOld + 1);
% StevensAgeOld = 27.778e7, equal to the 27.8 reported in Table 2  
Steven_kappaFCsOld_oneSigmaAbs = Steven(6)*1e6 * ...
                                (lambdaTotal_Steven*Steven_kappaFCsOld + lambdaEC_Steven);


% updating Steven et al. (1967) K-Ar data with new K IC and atomic weight data
wtPctK2O_to_mol40KperGramNew = (1/100)*(1/(2*39.0982+15.9994))*(2/1)*(1.1672e-4);
Steven_mol40KperGramNew = Steven(1)*wtPctK2O_to_mol40KperGramNew;
Steven_kappaFCsNew = Steven(3)/Steven_mol40KperGramNew;
% difference between kappaFCsOld and kappaFCsNew is 2.3%!  All from 40K/K value.
 
% back-calculate uncdrtainty in Steven_kappaFCsOld, take it as same uncertainty now.
% Steven explicity states (top of first column, page D54) no decay const. unct propagated.
% sigma_kappa = sigma_t * (lambdatotal*kappa + lambdaEC) from linear unct prop equation
Steven_kappaFCsNew_oneSigmaAbs = Steven(6)*1e6 * ...
                                (lambdaTotal_Steven*Steven_kappaFCsNew + lambdaEC_Steven);


%% Second, recreate Hurford and Hammerschmidt (1985)

% values from Garner et al. (1975) for SRM 985, from Table 11
r40K_K_Garner = 0.011672/100;
r40K_K_Garner_oneSigmaAbs = 0.000041/2/100;
r40K_K_Garner_oneSigmaRel = r40K_K_Garner_oneSigmaAbs/r40K_K_Garner;

wtPctK_to_mol40KperGram = (1/100)*(1/39.0983)*(r40K_K_Garner/1);
% (g K/100 g sample)*(mol K/g K)*(mol 40K/mol K)
% note: atomic mass of K from 1985 atomic masses

% n = PV/RT. @ STP, T = 273.15 K, P = 1 atm (pre 1982),
R = 82.057338; % cm^3 atm K^-1 mol^-1, used pre-1982, note 6 ppm chance since 1986.
%R = 83.144598; % cm^3 bar K^-1 mol^-1, used post-1982
cm3atSTP_to_molRad40ArPerGram = (1/R)*(1/273.15);

HandH_mol40KperGram = HandH(1)*wtPctK_to_mol40KperGram;
HandH_molRad40ArPerGram = HandH(3)*1e-6*cm3atSTP_to_molRad40ArPerGram;
HandH_kappaFCs = HandH_molRad40ArPerGram/HandH_mol40KperGram;
lambdaBeta_SandJ = 4.962e-10; % /year Steiger and Jager, 1977
lambdaEC_SandJ   = 0.581e-10; % /year Steiger and Jager, 1977
lambdaTotal_SandJ = lambdaBeta_SandJ + lambdaEC_SandJ;
HandHAge = 1/lambdaTotal_SandJ * log(lambdaTotal_SandJ/lambdaEC_SandJ*HandH_kappaFCs+1);
HandH_kappaFCs_oneSigmaAbs = HandH(6)*1e6 * ...
                                (lambdaTotal_SandJ*HandH_kappaFCs + lambdaEC_SandJ);
% difference between HandHAge/kappa and reported age is 0.1%, plausibly attibutable
% to roundoff errors in reporting.  Using age to calculate kappa for max precision
% and to match Renne et al 2010/11
HandH_reportedAge = (27.93+28.04)/2 * 10^6;
HandH_kappaFCs = lambdaEC_SandJ/lambdaTotal_SandJ * ...
                                           (exp(lambdaTotal_SandJ*HandH_reportedAge) - 1);

                                       
%% Third, calculate Jourdan and Renne FCS vs. other 'first principles' standards

% Hb3gr
% with no access to Turner et al. (1971), I'm stuck assuming data from Renne et al. 2010/11
R = 83.144598; % cm^3 bar K^-1 mol^-1, used post-1982
% note: using 39.0983 for mol/g K
kappaHb3gr =  JandR(1,3)*cm3atSTP_to_molRad40ArPerGram / ...
                 (JandR(1,1)*wtPctK_to_mol40KperGram);

% propagate uncertainties in kappaHb3gr, starting with just K and Ar determinations:
kappaHb3gr_oneSigmaAbsMeas = sqrt( ( cm3atSTP_to_molRad40ArPerGram / ...
                 (JandR(1,1)*wtPctK_to_mol40KperGram) )^2 * JandR(1,4)^2 + ...
                 (-JandR(1,3)*cm3atSTP_to_molRad40ArPerGram / ...
                 (JandR(1,1)^2*wtPctK_to_mol40KperGram))^2 * JandR(1,2)^2 );
% adding in the Garner (1975) reported uncertainty in 40K/K
kappaHb3gr_oneSigmaAbsW40K = sqrt( kappaHb3gr_oneSigmaAbsMeas^2 + ...
                                    (kappaHb3gr*r40K_K_Garner_oneSigmaRel)^2 );

% Reproduce the calculated K-Ar age and uncertainty of Hb3gr
Hb3grAge = 1/lambdaTotal_SandJ * log(lambdaTotal_SandJ/lambdaEC_SandJ*kappaHb3gr+1);
Hb3grAge_oneSigmaAbsMeas = sqrt( (1/(lambdaEC_SandJ+lambdaTotal_SandJ*kappaHb3gr))^2 * ...
                                  kappaHb3gr_oneSigmaAbsMeas^2 );
Hb3grAge_oneSigmaAbsW40K = sqrt( (1/(lambdaEC_SandJ+lambdaTotal_SandJ*kappaHb3gr))^2 * ...
                                  kappaHb3gr_oneSigmaAbsW40K^2 );

% relate back to FCs via the R reported in Jourdan and Renne (2007)
%kappaHb3gr_oneSigmaAbsMeas = 7e6 * ... % 1s abs for Hb3gr age from Jourdan and Renne 
%                                (lambdaTotal_SandJ*kappaHb3gr + lambdaEC_SandJ);
% note: this uncertainty of 1.4e-5 does not match any reported in Renne et al. 2010.
kappaFCs_Hb3gr = JandR(1,5)*kappaHb3gr; % kappa_FCs = kappa_Hb3gr * R_FCs/Hb3gr
kappaFCs_HB3gr_oneSigmaAbs = sqrt(kappaHb3gr^2*JandR(1,6)^2 + ...
                                               JandR(1,5)^2*kappaHb3gr_oneSigmaAbsMeas^2);

% GA-1550


%% compile re-calculated data
% data column 1: kappaFCs, measured or calibrated
% data column 2: one sigma absolute uncertainties

data(1,1) = Steven_kappaFCsOld; data(1,2) = Steven_kappaFCsOld_oneSigmaAbs;
data(2,1) = HandH_kappaFCs; data(2,2) = HandH_kappaFCs_oneSigmaAbs;

% I can't perfectly reproduce these calculations from the data given, so using values
% in Renne et al. 2010
data(3,:) = [1.665 0.023]*10^-3; % NL-25
data(4,:) = [1.641 0.010]*10^-3; % HB3gr    (I get [1.638 0.010])
data(5,:) = [1.639 0.008]*10^-3; % GHC-305
data(6,:) = [1.641 0.009]*10^-3; % GA-1550

kappaFCs = sum(data(:,1)./data(:,2).^2)/sum(1./data(:,2).^2);
kappaFCs_oneSigmaAbs = sqrt(1/sum(1./data(:,2).^2));
MSWD = 1/(size(data,1)-1)*sum( (data(:,1)-kappaFCs).^2./data(:,2).^2 );
dKappaFCs_dr40K_K = -kappaFCs/r40K_K_Garner;

% output vector:
kappaFCs = [kappaFCs, kappaFCs_oneSigmaAbs, dKappaFCs_dr40K_K];

end