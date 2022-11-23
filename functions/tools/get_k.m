%% get_k
% Check whether a new value of rate constant is provided in parameter
% Otherwise, use default value
%
%% Syntax
% |k = get_k(param,fieldName,default)|
%
%
%% Description
% |k = get_k(param,fieldName,default)| Description
%
%
%% Input arguments
% |param| -_STRUCTURE_- Data structure containing the new rate constants. |k=param.{fieldName}|
%
% |fieldName| -_STRING_- Name of the filed in |param| containing the new value of the constant
%
% |default| -_SCALAR_- If the field is absent, use |default as value for therate constant|
%
%
%% Output arguments
%
% |k| - _SCALAR_ - The rate constant
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function k = get_k(param,fieldName,default)

  if(isfield(param,fieldName))
    k = getfield(param,fieldName);
  else
    k=default;
  end
end
