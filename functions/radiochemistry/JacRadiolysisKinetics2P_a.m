%% JacRadiolysisKinetics
% Wrapper function to compute the Jacobian of the system of ODE defined in radiolysisKinetics.m
% Use this function as parameter 'Jacobian' in defining the parameters of the function odeset
%
%% Syntax
% |Jac = JacRadiolysisKinetics(t,C,param,TissueParam,biologyFLAG)|
%
%
%% Description
% |Jac = JacRadiolysisKinetics(t,C,param,TissueParam,biologyFLAG)| Description
%
%
%% Input arguments
% See radiolysisKinetics.m
%
%
%% Output arguments
%
% |Jac| -_SCALAR MATRIX_-  Jacobian matrix Jac(sp,j) = d_dydt(sp) / d_C(j) at time t, i.e. derivative of the reaction rate of species sp with respect to species j
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Jac = JacRadiolysisKinetics2P_a(t,C,param,TissueParam)

  [dydt , labels , Jac]= radiolysisKinetics2P_a(t,C,param,TissueParam);

end
