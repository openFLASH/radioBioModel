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


%on/off biological reactions
indexLegend = 1;
legendSTR = {};
legendSTR2 = {};

%Loop
for  NbPulses_i = 1:length(NbPulses)

      for doseRateIndex = 1:length(AvDoseRate)

          Period    = TotalDose./(NbPulses(NbPulses_i) .* AvDoseRate(doseRateIndex));  %Period (s)

          fprintf('Period : %f \n',Period)


          if (Period > PulseWidth(NbPulses_i))
            %Make computation only if the period is longer than the pulse duration

                [LOOHf(doseRateIndex, NbPulses_i) , AvDR  ] = getLOOHf(TotalDose , Period , PulseWidth(NbPulses_i) , NbPulses(NbPulses_i) , O2);

                legendSTR{end+1} =  ['<dD/dt> = ' num2str(AvDR,'%2.1g') ' Gy/s -- T = ' num2str(Period .* 1e3,'%2.1g') ' ms'];

                fprintf('Average dose rate = %f Gy/s \n',AvDR)
                fprintf('Dose %g Gy \n',TotalDose);
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
figure(201)
contour(DPP , AvDoseRate   , LOOHf , '-b', 'ShowText','on' )
ylabel('DR_a (Gy/s)')
xlabel('D_{pulse} (Gy)')
set(gca, 'YScale', 'log');
grid minor
title(['[LOOH]_f'])

AvDoseRate
idx = min(find(AvDoseRate > 100))
figure(202)
plot(DPP , squeeze(LOOHf(idx , :)),'o-')
xlabel('D_{pulse} (Gy)')
ylabel('[LOOH]_f (n-mol/l)')
title(['[LOOH]_f at DR_a = ' num2str(AvDoseRate(idx)) 'Gy/s'])
grid on


AvDoseRate
idx = max(find(DPP >= 4))
figure(203)
plot(AvDoseRate , squeeze(LOOHf(: , idx)),'o-')
xlabel('DR_a (Gy/s)')
ylabel('[LOOH]_f (n-mol/l)')
title(['[LOOH]_f at DPP = ' num2str(DPP(idx),'%3.0f') ' Gy'])
grid on
