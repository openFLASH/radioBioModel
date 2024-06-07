%% IsBeamOn
% Determine whether the time point is in BEAM ON or BEAM OFF period
%
%% Syntax
% |beam = IsBeamOn(t, param)|
%
%
%% Description
% |beam = IsBeamOn(t, param)| Description
%
%
%% Input arguments
% |t| - _SCALAR VECTOR_ - Time point at which the test is to be run
%
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants
% * |param.td| -_SCALAR_- (s) Total duration of the sequence of BEAM ON pulses
% * |param.t_on| -_SCALAR_- (s) For pulsed beam, duration of a single pulse
% * |param.T| -_SCALAR_- (s) Period separating 2 pulses
%
%% Output arguments
%
% |beam| - _BOOLEAN VECTOR_ - |beam()i = true| if the bema is On at time |t(i)|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function beam = IsBeamOn(t, param)
  beam = (mod(t, param.T) <= param.t_on) .* (t <= param.td);
end
