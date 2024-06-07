% Compute the concentration of different radical species as a function of time
% Vary the average dose rate and the dose per pulse
% Use the radio-kinetic model for the 2 phases model
%
%% Syntax
%
%REFERENCE
%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

clear
close all

folder = 'D:\programs\openREGGUI\REGGUI_userdata\radiokinetics'

store_y = {};
legendSTR = [];


colors = {'k','b','g','r','c','m','y'};
symbols = {'o', '*', '+',  '.', 'x', '_', '|', 'square'	, 'diamond'	, '^'	, 'v'	, '>'	, '<'	, 'pentagram', 'hexagram'};

%TotalDose =  10  ;% Gy
TotalDose =  14  ;% Gy
O2 = 50; %u-mol/l

%Beam parameter chosen to match experiment:
%Kevin Liu The FLASH effect is dependent on the dose per pulse or the mean dose rate for abdominal irradiation ESTRO 2024. https://www.biorxiv.org/content/10.1101/2024.04.19.590158v1.full.pdf
PeakDR = 1.7e6 %Gy/s Peak dose rate
NbPulses = 1:8 %Make sure that there is an integral number of complete pulses
                %This parameter controls the dose per pulse

DPP = TotalDose ./ NbPulses %Gy Dose per pulse
PulseWidth = DPP ./ PeakDR; %s In the Mobetron, the dose per pulse is controlled by increasing the width of the pulse and jkeeping the peak DR constant

idx = 1:10;
AvDoseRate = 2.^idx %Gy/s between 2 and 1020 Gy/s

%Save the simulation parameters
data.TotalDose = TotalDose;
data.AvDoseRate = AvDoseRate ;
data.PulseWidth = PulseWidth ;
data.NbPulses = NbPulses;
data.O2 = O2;


%Rate constants
%------------------
[kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants();
param = getDefaultParam();
kName = {  'kbr2' , 'kb3' , 'kLOOself' , 'kb8' , 'kbr'};
kValue = [kbr2 , kb3,  kLOOself , kb8, kbr];
param = set_k_in_param(param,kValue,kName); %Update the value of the rate constant


%on/off biological reactions
indexLegend = 1;
legendSTR = {};
legendSTR2 = {};

%Loop
for  NbPulses_i = 1:length(NbPulses)
      for doseRateIndex = 1:length(AvDoseRate)

          param.R0   = AvDoseRate(doseRateIndex); %Average dose rate
          param.T    = TotalDose./(NbPulses(NbPulses_i) .* AvDoseRate(doseRateIndex));  %Period (s)
          param.t_on = PulseWidth(NbPulses_i); %s Duration of a single pulse

          fprintf('Period : %f \n',param.T)

          if (param.T > param.t_on)
            %Make computation only if the period is longer than the pulse duration

                legendSTR{end+1} =  ['<dD/dt> = ' num2str(param.R0,'%2.1g') ' Gy/s -- T = ' num2str(param.T .* 1e3,'%2.1g') ' ms'];

                fprintf('Computing for average dose rate %1.3g Gy/s .... \n',param.R0);
                fprintf('Period  : %g s \n',param.T)
                fprintf('Beam ON : %1.3g s \n',param.t_on)

                param.td = TotalDose ./ param.R0; % s Beam ON time. Set the beam time to deliver the same total dose for all dose rate
                fprintf('Computing for dose %g Gy => delivery time %g s \n',TotalDose,param.td);
                fprintf('Number of pulses : %d \n',ceil(param.td ./ param.T))

                %Initial concentrations
                %=======================
                [dydt , labels]= radiolysisKinetics2P_a();
                Nbconcentrations = length(labels);
                fprintf('Number of tracked concentrations: %d \n',Nbconcentrations);
                labels

                [y0 , t0] = getY0(param , @radiolysisKinetics2P_a); %Radical concentration at begining of homogeneous chemical phase
                y0(2)  = O2; %u-mol/l Oxygen inital concentration

                if (strcmp(func2str(param.R) , 'pulsedBeam'))
                  fprintf('Pulsed beam -- Pulse width = %g s \n',param.t_on);
                  opts = odeset('NonNegative',1:length(y0),'InitialStep',param.t_on.*1e-4,'Jacobian',@JacRadiolysisKinetics2P_a);
                else
                  fprintf('Not a pulsed beam -- Default ode time step \n')
                  opts = odeset('NonNegative',1:length(y0),'Jacobian',@JacRadiolysisKinetics2P_a);
                end
                Tstart = datetime;
                fprintf('Computation starts at %s \n',Tstart);
                [t,y] = pulsedODE(@odefcnComprehensive, t0 , y0,opts,param,param);
                Tend = datetime;
                fprintf('Computation ends at %s \n',Tend);
                fprintf('Duration : %s \n',Tend-Tstart);

                % Store all computation results per O2/totalDose/doseRate
                store_y{doseRateIndex}{NbPulses_i} = [t,y];

                [~ , Rp] = param.R(0,param);
                fprintf('Peak dose rate = %2.3g Gy/s \n',Rp);

                %Display graph for all species
                %=============================
                %displayGraphSpecies(t, y , doseRateIndex , legendSTR , kb8_i);

                %Get final LOOH concentration
                %=====================================
                LOOHf(doseRateIndex, NbPulses_i) =  y(end,10); %uM
                fprintf('Average dose rate = %f Gy/s \n',param.R0)
                fprintf('Dose %g Gy => delivery time %g s \n',TotalDose,param.td);
                fprintf('[LOOH]f = %g uM \n',LOOHf(doseRateIndex, NbPulses_i))
                fprintf('DONE \n')
          else
            %No computation done
            fprintf('SKIP \n')
            LOOHf(doseRateIndex, NbPulses_i) =  NaN;

          end %if (param.T > param.t_on)
      end %for doseRateIndex

      legendSTR2{end+1} =  ['D_{pulse} = ' num2str(TotalDose ./ NbPulses(NbPulses_i), '%2.2g') ' Gy'];

      figure(200)
      semilogx(AvDoseRate,squeeze(LOOHf(:, NbPulses_i)), ['-' symbols{mod(NbPulses_i,length(symbols))+1} colors{mod(NbPulses_i,length(colors))+1}]);
      hold on
      xlabel('Average dose rate (Gy/s)')
      ylabel('[LOOH]_f (\mu M)')
      title(['[LOOH]_f'  ])
      legend(legendSTR2 , 'Location' , 'eastoutside')
      grid on
      drawnow


end %for NbPulses_i


DPP = TotalDose ./ NbPulses;
figure(1)
contour(DPP , AvDoseRate   , LOOHf , '-b', 'ShowText','on' )
ylabel('DR_a (Gy/s)')
xlabel('D_{pulse} (Gy)')
set(gca, 'YScale', 'log');
grid minor
title(['[LOOH]_f'])

AvDoseRate
idx = min(find(AvDoseRate > 100))
figure(2)
plot(DPP , squeeze(LOOHf(idx , :)),'o-')
xlabel('D_{pulse} (Gy)')
ylabel('[LOOH]_f (n-mol/l)')
title(['[LOOH]_f at DR_a = ' num2str(AvDoseRate(idx)) 'Gy/s'])
grid on


AvDoseRate
idx = max(find(DPP >= 4))
figure(3)
plot(AvDoseRate , squeeze(LOOHf(: , idx)),'o-')
xlabel('DR_a (Gy/s)')
ylabel('[LOOH]_f (n-mol/l)')
title(['[LOOH]_f at DPP = ' num2str(DPP(idx),'%3.0f') ' Gy'])
grid on

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
  fprintf('conversion factor kr-ge = %g \n',ge2kr(1));
  param.pH = 7; %pH of the extra vascular tissue. NB: assume buffered solution [H+] and [OH-] are constant

  param.R = @pulsedBeam;
  %param.R = @constantBeam;


end
