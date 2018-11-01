function S1X = KwonPlusRenneTerms(x, datePairs, deltaEruptionZircon, K40decay, kappaFCs)
%KWONETAL2002_OPTFUN Optimization function and approach from Kwon et al. 2002
%   Setup using equations 4 and 5 of Renne et al. 2010
%   datePairs are [t_zircon (Ma), 1s abs, R_FCs, 1s abs]
%   x = [lambdaEC, lambdaBeta, kappaFCs, R_i=1_to_n]

lambdaECStar = x(1);
lambdaBetaStar = x(2);
lambdaTotalStar = lambdaECStar + lambdaBetaStar;
kappaFCsStar = x(3); % don't worry about tossing out other values here
Rstar = x(4:end)';

% equation 4
tauStar = 1/lambdaTotalStar * log(lambdaTotalStar/lambdaECStar * kappaFCsStar*Rstar + 1); 

t  = datePairs(:,1) * 1e6; % convert from Ma to a for decay constants in a^-1
st = datePairs(:,2) * 1e6; 
R  = datePairs(:,3);
sR = datePairs(:,4);

t = t - deltaEruptionZircon(1);
st = sqrt(st.^2 + deltaEruptionZircon(2)^2);

SX = sum( (1/2)*((tauStar-t)./st).^2 + (1/2)*((Rstar-R)./sR).^2 ); % equation 5

S1X = SX + ( (x(1)-K40decay(1)) / K40decay(2) )^2 + ...
           ( (x(2)-K40decay(3)) / K40decay(4) )^2 + ...
           ( (x(3)-kappaFCs(1)) / kappaFCs(2) )^2;             % equation 6

end

