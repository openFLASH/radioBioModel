%% modelAUCvsX_DR
% Compute the area under the curve (AUC) of the alkyl-peroxyl radicals (ROO) for a given dose rate and a second parameter (x).
%
%   AUCmodel = a(x)+b(x) ./ sqrt(doseRate);
%
% where the function a(x) and b(x) are interpolated from a table of value as a function of parameter (x).
% for example, x can be the dose or the oxygen concentration.
%
%% Syntax
% |AUCmodel = modelAUCvsX_DR(param,X,doseRate)|
%
%
%% Description
% |AUCmodel = modelAUCvsX_DR(param,X,doseRate)| compute AUC
%
%
%% Input arguments
% |param| - _SCALAR MATRIX_ - The table giving the value of the coefficients a =|param(i,1)| (in uM.s) and b =|param(i,2)| (in uM.s.Gy^0.5) at parameter X =|param(i,3)|
%
% |X| - _SCALAR VECTOR_ - The AUC(i) (in uM.s) is computed for parameter |X(i)| and dose rate |doseRate(i)| (Gy/s)
%
% |doseRate| - _SCALAR VECTOR_ - The AUC(i) is computed for parameter |X(i)| and dose rate |doseRate(i)|
%
%
%% Output arguments
%
% |AUCmodel|  - _SCALAR VECTOR_ - The AUC(i) (in uM.s) is computed for dose |dose(i)| and dose rate |doseRate(i)|
%
% |a| - _SCALAR VECTOR_ - The coefficients a =|param(i,1)| (in uM.s) at parameter X =|X(i)| (in Gy)
%
% |b| - _SCALAR VECTOR_ - The coefficients b =|param(i,2)| (in uM.s.Gy^0.5) at parameter X =|X(i)| (in Gy)
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [AUCmodel , a , b ]= modelAUCvsX_DR(param,X,doseRate)
  aM = param(:,1);
  bM = param(:,2);
  Dref = param(:,3);

  %Interpolate the a,b parameter for the required parameter X
  pa = polyfit(Dref,squeeze(aM),4);
  pb = polyfit(Dref,squeeze(bM),4);

  a = polyval(pa,X);
  b = polyval(pb,X);

  AUCmodel = a + b ./ sqrt(doseRate);
end
