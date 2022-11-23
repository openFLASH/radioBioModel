%% constantBeam
% Determine whether the time point is in BEAM ON or BEAM OFF for a constant dose rate
%
%% Syntax
% |[doseRate , Rp] = constantBeam(t,param)|
%
%
%% Description
% |[doseRate , Rp] = constantBeam(t,param)| Description
%
%
%% Input arguments
% |t| - _SCALAR VECTOR_ - Time point at which the test is to be run
%
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants
% * |param.td| -_SCALAR_- (s) Total duration of the sequence of BEAM ON pulses
% * |param.R0| -_SCALAR_- (Gy/s) Dose rate
%
%% Output arguments
%
% |doseRate| - _SCALAR VECTOR_ - |doseRate(i)| is the dose rate (Gy/s) at time |t(i)|
%
% |Rp| - _SCALAR_ - Average dose rate (Gy/s)
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [doseRate , Rp] = constantBeam(t,param)

Rp = param.R0; %Gy/s;
if (t <= param.td) %Beam on time
  doseRate = param.R0; %Gy/s
else
  doseRate = 0;
end

end
