%% set_k_in_param
% Update value of rate constants in structure |param|
%
%% Syntax
% |param = set_k_in_param(param,k,kName)|
%
%
%% Description
% |param = set_k_in_param(param,k,kName)| Add to |param| the fields defined in |kanme| with the value defined in |k|
%
%
%% Input arguments
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants. See |radiolysisKinetics|
%
% |k| - _SCALAR VECTOR_ - |k| Value of the i-th field to add to param
%
% |kName| - _CELL VECTOR_ - |kName(i)| Name of the i-th field to add to |param|
%
%% Output arguments
%
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants. Updated structure
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function param = set_k_in_param(param,k,kName)
  for ki=1:length(k)
    param = setfield(param,kName{ki},k(ki));
    fprintf('%s = %3.3g :: ',kName{ki},k(ki))
  end
  fprintf('\n')
end
