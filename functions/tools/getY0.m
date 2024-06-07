%% getY0
% Compute the concentration at the start of the simulation.
% The simulation starts at the begining of the homogeneous chemical phase (1us post irradiaiton)
%
%% Syntax
% |[y0 , t0] = getY0(param , @radioModel)|
%
%
%% Description
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants
%
% |@radioModel| -_FUNCTION POINTER_- POinter to the radio kinetic model function
%
%% Input arguments
% |im1| - _STRING_ -  Name
%
%
%% Output arguments
%
% |y0| - _SCALAR VECTOR_ - Inital concentration of the chemica lspecies of the model at 1us post irradiation
%
% |t0| -_SCALAR_- Time (s) at which the |y0| are computed
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [y0 , t0] = getY0(param , radioModel)

  %[~ , labels]= radiolysisKinetics2P_a();
  [~ , labels]= radioModel();

  G = zeros(numel(labels) , 1); % mol/l/Gy
  G(find(strcmp(labels, 'e^-_{aq}' ))) = param.kre;
  G(find(strcmp(labels, 'OH^.' ))) = param.krOHr;
  G(find(strcmp(labels, 'H^.' ))) = param.krHr;
  G(find(strcmp(labels, 'H_2' ))) = param.krH2;
  G(find(strcmp(labels, 'H_2O_2' ))) = param.krH2O2;

  t0 = min(1e-6 , param.t_on); %s Integrate delivered dose either to the end of pulse (if less than 1us) OR to the start of the homogeneous chemistry phase
  D = integral(@(x) param.R(x,param) , 0 , t0); %Total dose Gy delivered up to time t0
  t0 = 1e-6;  %s The simulation starts at the homogeneous chemistry phase
  y0 = G .* D; %mol/l of the different chemicals at the end of the

end
