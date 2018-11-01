function K40decay = Kdecay()

%% compile 40K decay constant data used by Min et al (2000)

atomicWt_K = 39.0983; % current as of 2016 IUPAC (Meija et al., 2016)
atomicWt_K_oneSigmaAbs = 0.00012/2; % standardized units used for standard atomic weights
ratio_40K_K =  1.17e-4; % ENSDF, Audi et al., 1997
ratio_40K_K_oneSigmaAbs = (0.02e-4)/2;

avogadrosNumber = 6.0221367e23; % mol^-1  
avogadrosNumber_oneSigmaAbs = (0.0000072e23)/2; % mol^-1
% Current estimate Oct 2018: 6.022140857e23, difference = -0.7 ppm

numberOfSecondsInOneTropicalYear =  3.155693e7;
% Jaffey uses 5.2595e5 minutes = 3.1557e7 seconds, difference = -2.2 ppm

% Activity averages and uncertainties for electron capture and beta decay
% reported in Table 1 of Min et al. (2000)
activityEC_EandVdL = 3.31; % g^-1 sec^-1
activityEC_EandVdL_oneSigmaAbs = 0.06/2; % g^-1 sec^-1
activityBeta_EandVdL = 27.89; % g^-1 sec^-1
activityBeta_EandVdL_oneSigmaAbs = 0.30/2; % g^-1 sec^-1


%% Calculate decay constants from activity data
% from equation 1 of Min et al. (2000), lambda in units of year^-1

activityToLambda = atomicWt_K * numberOfSecondsInOneTropicalYear / ...
                   (ratio_40K_K * avogadrosNumber); % W*Y/(f*N_A)               
lambdaEC   = activityEC_EandVdL   * activityToLambda; % lambda = A * W*Y/(f*N_A)
lambdaBeta = activityBeta_EandVdL * activityToLambda;


%% Uncertainty propagation

% electron capture
dLambdaEC_dActivity = activityToLambda; % W*Y/(f*N_A)
dLambdaEC_dAtomicWt = activityEC_EandVdL * numberOfSecondsInOneTropicalYear / ...
                      (ratio_40K_K * avogadrosNumber); % A*Y/(f*N_A)
dLambdaEC_d40K_K = -activityEC_EandVdL*atomicWt_K*numberOfSecondsInOneTropicalYear /...
                    (ratio_40K_K^2 * avogadrosNumber); % -A*W*Y/(f^2*N_A)
dLambdaEC_dNA = -activityEC_EandVdL*atomicWt_K*numberOfSecondsInOneTropicalYear /...
                    (ratio_40K_K * avogadrosNumber^2); % -A*W*Y/(f*N_A^2)

lambdaEC_variance = dLambdaEC_dActivity^2 * activityEC_EandVdL_oneSigmaAbs^2 + ...
                    dLambdaEC_dAtomicWt^2 * atomicWt_K_oneSigmaAbs^2 + ...
                    dLambdaEC_d40K_K^2    * ratio_40K_K_oneSigmaAbs^2 + ...
                    dLambdaEC_dNA^2       * avogadrosNumber_oneSigmaAbs^2;
lambdaEC_oneSigmaAbs = sqrt(lambdaEC_variance);

% beta- decay
dLambdaBeta_dActivity = activityToLambda; % W*Y/(f*N_A)
dLambdaBeta_dAtomicWt = activityBeta_EandVdL * numberOfSecondsInOneTropicalYear / ...
                      (ratio_40K_K * avogadrosNumber); % A*Y/(f*N_A)
dLambdaBeta_d40K_K = -activityBeta_EandVdL*atomicWt_K*numberOfSecondsInOneTropicalYear/...
                    (ratio_40K_K^2 * avogadrosNumber); % -A*W*Y/(f^2*N_A)
dLambdaBeta_dNA = -activityBeta_EandVdL*atomicWt_K*numberOfSecondsInOneTropicalYear /...
                    (ratio_40K_K * avogadrosNumber^2); % -A*W*Y/(f*N_A^2)

lambdaBeta_variance = dLambdaBeta_dActivity^2 * activityBeta_EandVdL_oneSigmaAbs^2 + ...
                    dLambdaBeta_dAtomicWt^2 * atomicWt_K_oneSigmaAbs^2 + ...
                    dLambdaBeta_d40K_K^2    * ratio_40K_K_oneSigmaAbs^2 + ...
                    dLambdaBeta_dNA^2       * avogadrosNumber_oneSigmaAbs^2;
lambdaBeta_oneSigmaAbs = sqrt(lambdaBeta_variance);

covEC_Beta = dLambdaEC_dAtomicWt*dLambdaBeta_dAtomicWt*atomicWt_K_oneSigmaAbs^2 + ...
             dLambdaEC_d40K_K   *dLambdaBeta_d40K_K   *ratio_40K_K_oneSigmaAbs^2 + ...
             dLambdaEC_dNA      *dLambdaBeta_dNA      *avogadrosNumber_oneSigmaAbs^2;
rhoEC_Beta = covEC_Beta/(lambdaEC_oneSigmaAbs*lambdaBeta_oneSigmaAbs);

lambdaTotal = lambdaEC + lambdaBeta;
lambdaTotal_variance = lambdaEC_variance + lambdaBeta_variance + 2*covEC_Beta;
lambdaTotal_oneSigmaAbs = sqrt(lambdaTotal_variance);


%% Output
% note: revision in Renne et al. (2011) comment discards liquid scintillation data.

K40decay = [lambdaEC   lambdaEC_oneSigmaAbs ...
            lambdaBeta lambdaBeta_oneSigmaAbs ...
            rhoEC_Beta ...
            lambdaTotal lambdaTotal_oneSigmaAbs];