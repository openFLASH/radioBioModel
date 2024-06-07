%% displayGraphSpecies
% Compute the [LOOH]f for the different specified experimental conditions
% All input vector have the smae length. Each element of the vectors define one experimental condition.
%
%% Syntax
% |[LOOHf , AvDoseRate , PkDoseRate ] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue)|
%
% |[LOOHf , AvDoseRate , PkDoseRate ] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 )|
%
%
%% Description
% |[LOOHf , AvDoseRate , PkDoseRate ] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue)| Description
%
%
%% Input arguments
% |TotalDose| - _SCALAR VECTOR_ - (Gy) Total delivered dose
%
% |Period| - _SCALAR VECTOR_ -  (s) Period of the pulsed beam
%
% |PulseWidth| - _SCALAR VECTOR_ - (s) Duration of each pulse
%
% |NbPulses| - _SCALAR VECTOR_ - Number of pulses
%
% |O2| - _SCALAR VECTOR_ - (u-mol/l) Oxygen concentration
%
% |kValue| - _SCALAR VECTOR_ - [OPTIONAL] [kbr2 , kb3,  kLOOself , kb8 , kbr] Value of rate constants to use in the model
%
% |verbose| -_BOOLEAN_- [OPTIONAL. Default: true]
%
%% Output arguments
%
% |LOOHf| - _SCALAR VECTOR_ - (u-mol/l) Concentration of LOOH at the end of the experiment
%
% |AvDoseRate| - _SCALAR VECTOR_ - (Gy/s) Average dose rate of experiment
%
% |PkDoseRate| - _SCALAR VECTOR_ - (Gy/s) Peak dose rate of experiment
%
% |O2f| - _SCALAR VECTOR_ - (u-mol/l) (u-mol/l) Concentration of LOOH at the end of the experiment
%
%% Contributors

function [LOOHf , AvDoseRate , PkDoseRate , O2f] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue , verbose)

  if nargin < 7
    verbose = true;
  end

  %Rate constants
  %------------------
  if nargin < 6
    %Use default rate constatns
    [kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants();
    kValue = [kbr2 , kb3,  kLOOself , kb8 , kbr];
  elseif isempty(kValue)
    [kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants();
    kValue = [kbr2 , kb3,  kLOOself , kb8 , kbr];
  else
    %kValue is provided as input
  end

  %Manually change the value of some rate constants in the model
  param = getDefaultParam();
  kName = {  'kbr2' , 'kb3' , 'kLOOself' , 'kb8' , 'kbr'};
  param = set_k_in_param(param,kValue,kName); %Update the value of the rate constant

  %Compute average dose ratez
  AvDoseRate =  TotalDose ./ (NbPulses .* Period); %Gy/s %average dose rate

  legendSTR = {};

  %Loop for every simulation
  for idx = 1:numel(O2)

      param.R0   = AvDoseRate(idx); %Average dose rate
      param.T    = Period(idx); %s
      param.t_on = PulseWidth(idx); %s Duration of a single pulse
      param.td = TotalDose(idx) ./ param.R0; % s Beam ON time. Set the beam time to deliver the same total dose for all dose rate
      [~ , Rp] = param.R(0,param);
      PkDoseRate(idx) = Rp;

      if verbose
        legendSTR{end+1} =  ['<dD/dt> = ' num2str(param.R0,'%2.1g') ' Gy/s -- (dD/dt)_P = ' num2str(Rp,'%2.1g') ' Gy/s'];
        fprintf('Computing for average dose rate %1.3g Gy/s .... \n',param.R0);
        fprintf('Period  : %g s \n',param.T)
        fprintf('Beam ON : %1.3g s \n',param.t_on)
        fprintf('Computing for dose %g Gy => delivery time %g s \n',TotalDose(idx),param.td);
        fprintf('Number of pulses : %d \n',ceil(param.td ./ param.T))
        fprintf('Peak dose rate = %2.3g Gy/s \n',Rp);
      end


      %Initial concentrations
      %=======================
      [dydt , labels]= radiolysisKinetics2P_a();
      Nbconcentrations = length(labels);
      if verbose
        fprintf('Number of tracked concentrations: %d \n',Nbconcentrations);
        labels
      end

      [y0 , t0] = getY0(param , @radiolysisKinetics2P_a); %Radical concentration at begining of homogeneous chemical phase
      y0(2)  = O2(idx); %u-mol/l Oxygen inital concentration

      if (strcmp(func2str(param.R) , 'pulsedBeam'))
        if verbose
          fprintf('Pulsed beam -- Pulse width = %g s \n',param.t_on);
        end
        opts = odeset('NonNegative',1:length(y0),'InitialStep',param.t_on.*1e-4,'Jacobian',@JacRadiolysisKinetics2P_a);
      else
        if verbose
          fprintf('Not a pulsed beam -- Default ode time step \n')
        end
        opts = odeset('NonNegative',1:length(y0),'Jacobian',@JacRadiolysisKinetics2P_a);
      end
      Tstart = datetime;
      if verbose
        fprintf('Computation starts at %s \n',Tstart);
      end
      [t,y] = pulsedODE(@odefcnComprehensive, t0 , y0,opts,param,param);
      Tend = datetime;
      if verbose
        fprintf('Computation ends at %s \n',Tend);
        fprintf('Duration : %s \n',Tend-Tstart);

        %Display graph for all species
        %=============================
        displayGraphSpecies(t, y , idx , legendSTR , idx);
      end

      %Get final LOOH concentration
      %=====================================
      [~ , labels]= radiolysisKinetics2P_a();
      LOOHidx = find(strcmp(labels, 'LOOH')); %find the index for [O2]
      LOOHf(idx) =  y(end,LOOHidx); %uM
      if verbose
        fprintf('Average dose rate = %f Gy/s \n',param.R0)
        fprintf('Dose %g Gy => delivery time %g s \n',TotalDose(idx),param.td);
        fprintf('[LOOH]f = %g uM \n',LOOHf(idx))
      end


      %Get final O2 concentration
      %============================
      [~ , labels]= radiolysisKinetics2P_a();
      O2idx = find(strcmp(labels, 'O_2')); %find the index for [O2]
      O2f(idx) =  y(end,O2idx); %uM
      if verbose
        fprintf('[O2]f = %g uM \n',O2f(idx))
        fprintf('DONE \n')
      end

  end %for idx

end


%==================================
% system of Ordinary differential Equation
% Kinetics of oxygen depletion by aqueous electrons
%==================================

function dydt = odefcnComprehensive(t,y,param,TissueParam)
  dydt = radiolysisKinetics2P_a(t,y,param,TissueParam);
end


%----------------------------------
% Define the default value for the paraemeters
%----------------------------------
function param = getDefaultParam()
  param = paramRadiolytic();
  param.pH = 7; %pH of the extra vascular tissue. NB: assume buffered solution [H+] and [OH-] are constant
  param.R = @pulsedBeam;
  %param.R = @constantBeam;


end
