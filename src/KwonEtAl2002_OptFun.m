function SX = KwonEtAl2002_OptFun(x, datePairs, deltaEruptionZircon)
%KWONETAL2002_OPTFUN Optimization function and approach from Kwon et al. 2002
%   Setup using equations 4 and 5 of Renne et al. 2010
%   datePairs are [t_zircon (Ma), 1s abs, R_FCs, 1s abs]
%   x = [lambdaEC, lambdaBeta, kappaFCs, R_i=1_to_n]

lambdaEC = x(1);
lambdaBeta = x(2);
lambdaTotal = lambdaEC + lambdaBeta;
kappaFCs = x(3); % don't worry about tossing out other values here
Rstar = x(4:end)';
tau = 1/lambdaTotal * log(lambdaTotal/lambdaEC * kappaFCs*Rstar + 1); % equation 4

t  = datePairs(:,1) * 1e6; % convert from Ma to a for decay constants in a^-1
st = datePairs(:,2) * 1e6; 
R  = datePairs(:,3);
sR = datePairs(:,4);

t = t - deltaEruptionZircon(1);
st = sqrt(st.^2 + deltaEruptionZircon(2)^2);

SX = sum( (1/2)*((tau-t)./st).^2 + (1/2)*((Rstar-R)./sR).^2 ); % equation 5

end

