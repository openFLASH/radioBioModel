%% getConc
% Compute the concentration of the species for the different specified experimental conditions
%
%% Syntax
% |[LOOHf , AvDoseRate , PkDoseRate ] = getConc(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue)|
%
% |[LOOHf , AvDoseRate , PkDoseRate ] = getConc(TotalDose , Period , PulseWidth , NbPulses , O2 )|
%
%
%% Description
% |[LOOHf , AvDoseRate , PkDoseRate ] = getConc(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue)| Description
%
%
%% Input arguments
% |TotalDose| - _SCALAR_ - (Gy) Total delivered dose
%
% |Period| - _SCALAR_ -  (s) Period of the pulsed beam
%
% |PulseWidth| - _SCALAR_ - (s) Duration of each pulse
%
% |NbPulses| - _SCALAR_ - Number of pulses
%
% |O2| - _SCALAR_ - (u-mol/l) Oxygen concentration
%
% |kValue| - _SCALAR VECTOR_ - [OPTIONAL] [kbr2 , kb3,  kROOself , kb8 , kbr] Value of rate constants to use in the model
%
% |verbose| -_BOOLEAN_- [OPTIONAL. Default: true]
%
%% Output arguments
%
% |y0| - _SCALAR MATRIX_ - Concentration of the chemical species of the model at 1us post irradiation
%
% |t| -_SCALAR VECTOR_- Time (s) at which the |y| are computed
%
% |labels| -_CELL VECTOR_- Name of the tracked species
%
%% Contributors

function [t, y , labels] = getConc(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue , verbose)

  if nargin < 7
    verbose = true;
  end

  %Rate constants
  %------------------

  if nargin < 6
    %Use default rate constatns
    [kbr2 , kROOself , kb3 , kb8 , kbr] = getRateConstants();
    kValue = [kbr2 , kb3,  kROOself , kb8 , kbr];
  elseif isempty(kValue)
    [kbr2 , kROOself , kb3 , kb8 , kbr] = getRateConstants();
    kValue = [kbr2 , kb3,  kROOself , kb8 , kbr];
  else
    %kValue is provided as input
  end

  %Manually change the value of some rate constants in the model
  param = getDefaultParam();
  kName = {  'kbr2' , 'kb3' , 'kROOself' , 'kb8' , 'kbr'};
  param = set_k_in_param(param,kValue,kName); %Update the value of the rate constant

  %Compute average dose ratez
  AvDoseRate =  TotalDose ./ (NbPulses .* Period); %Gy/s %average dose rate

  legendSTR = {};

  %Simulation
      param.R0   = AvDoseRate; %Average dose rate
      param.T    = Period; %s
      param.t_on = PulseWidth; %s Duration of a single pulse
      param.td = TotalDose ./ param.R0; % s Beam ON time. Set the beam time to deliver the same total dose for all dose rate
      [~ , Rp] = param.R(0,param);
      PkDoseRate = Rp;

      if verbose
        legendSTR{end+1} =  ['<dD/dt> = ' num2str(param.R0,'%2.1g') ' Gy/s -- (dD/dt)_P = ' num2str(Rp,'%2.1g') ' Gy/s'];
        fprintf('Computing for average dose rate %1.3g Gy/s .... \n',param.R0);
        fprintf('Period  : %g s \n',param.T)
        fprintf('Beam ON : %1.3g s \n',param.t_on)
        fprintf('Computing for dose %g Gy => delivery time %g s \n',TotalDose,param.td);
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
      y0(2)  = O2; %u-mol/l Oxygen inital concentration

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
        displayGraphSpecies(t, y , 1 , legendSTR , 1);
      end


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
