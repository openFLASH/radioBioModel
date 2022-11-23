%% fitAUCROOr
% Fit the model of the area under the curve (AUC) of the alkyl-peroxyl radicals (ROO) as a function of  dose rate.
%   AUCmodel = a + b ./ sqrt(doseRate)
%
%% Syntax
% |param = fitAUCROOr(doseRate,INTRooR,symbol)|
%
%
%% Description
% |param = fitAUCROOr(doseRate,INTRooR,symbol)| Description
%
%
%% Input arguments
% |doseRate| - _SCALAR VECTOR_ - |doseRate| The i-th dose rate (Gy/s) at which the AUC is computed
%
% |INTRooR|  - _SCALAR VECTOR_ - |INTRooR| The AUC at the i-th dose rate
%
% |symbol| -_STRING_- [OPTIONAL] If present, a plot of the fit is shown. The symbol is used to represent the AUC values
%
% |colour| -_STRING_- [OPTIONAL] If present (and if |symbol is present|, define the colour of the symbol and line on the plot)
%
%% Output arguments
%
% |param| -_SCALAR VECTOR_- [a,b] parameter of the fitted curve
%
%
%% Contributors
% Authors : R.Labarbe (open.reggui@gmail.com)
function param = fitAUCROOr(doseRate,INTRooR,symbol,colour)


  %Fit the ROOr vs dose rate curve
  %===============================
  options = optimset('Display','iter');
  b = (INTRooR(1) - INTRooR(end))./(1./sqrt(doseRate(1)) - 1./sqrt(doseRate(end)));
  a= INTRooR(1)-b./sqrt(doseRate(1));
  param0 = [a,b]
  [param,fval,exitflag] = fminsearch(@objective,param0,options,INTRooR', doseRate);

  fprintf('a = %f \n',param(1))
  fprintf('b = %f \n',param(2))

  if(nargin > 2)
    DR = power(10,min(log10(doseRate)):0.1:max(log10(doseRate)));
    AUCmodel = modelAUCvsDR(param,DR)
    figure(20)
    if(nargin >3)
      symbol1 = [symbol,colour];
      symbol2=['-',colour]
    else
      symbol1 = symbol;
      symbol2='-';
    end
    semilogx(doseRate , INTRooR,symbol1)
    hold on
    semilogx(DR , AUCmodel,symbol2)
  end

end


%============================
% Objective function to minimize
%============================
function gof = objective(param, AUC , doseRate)
  AUCmodel = modelAUCvsDR(param,doseRate);
  gof = sum((AUC-AUCmodel).^2);
end
