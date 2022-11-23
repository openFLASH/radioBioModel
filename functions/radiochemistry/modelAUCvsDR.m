%% modelAUCvsX_DR
% Compute the area under the curve (AUC) of the alkyl-peroxyl radicals (ROO) for a given dose rate.
%
%   AUCmodel = a + b ./ sqrt(doseRate);
%
%
%% Syntax
% |AUCmodel = modelAUCvsDR(param,doseRate)|
%
%
%% Description
% |AUCmodel = modelAUCvsDR(param,doseRate)| compute AUC
%
%
%% Input arguments
% |param| - _SCALAR VECTOR_ - [a,b] The  coefficients a  (in uM.s) and b  (in uM.s.Gy^0.5)
%
% |doseRate| - _SCALAR VECTOR_ - The AUC(i) is computed for  dose rate |doseRate(i)|
%
%
%% Output arguments
%
% |AUCmodel|  - _SCALAR VECTOR_ - The AUC(i) (in uM.s) is computed for dose rate |doseRate(i)|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function AUCmodel = modelAUCvsDR(param,doseRate)
  a = param(:,1);
  b = param(:,2);
  AUCmodel = a+b./sqrt(doseRate);
end
