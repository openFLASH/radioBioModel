%==================================
% Driving pulsed dose rate function
%
% Input
% |t| -_SCALAR VECTOR_- Time point (s) at which the instantaneous dose rate (Gy/s) is computed
% |param| -_Structure_- Parameters describing the beam
%  * |param.R0| -_SCALAR_- Average dose rate (Gy/s)
%  * |param.t_on| -_SCALAR_- within 1 cycle, time during which the beam is on: |param.t_on| = |param.T| .* duty cycle
%  * |param.T| -_SCALAR_- Period (s) of one pulse sequence
%==================================
function [pulseDoseRate , Rp] = pulsedBeam(t, param)

    dutyCycle = param.t_on ./ param.T;
    Rp = param.R0 ./ dutyCycle; %Gy/s, convert from average dose rate into peak dose rate

    %pulseDoseRate = IsBeamOn(t, param) .* Rp;
    pulseDoseRate = (mod(t, param.T) <= param.t_on) .* (t <= param.td) .* Rp;

end
