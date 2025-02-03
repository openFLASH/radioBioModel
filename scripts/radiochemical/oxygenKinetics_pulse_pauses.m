% Compute the concentration of different radical species as a function of time
% Vary the 4 parameters of the pulse
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

NbDelivery =         [1   ,2  , 3   , 4  , 6] ; % The dose is splitted in several deliveries
TotalDose  = 39.9 ./ NbDelivery;% Gy
O2         =         [50 , 50 , 50 , 50 , 50  ]; %u-mol/l
NbPulses   =         [1  , 1  , 1  , 1  , 1 ];
PulseWidth = TotalDose ./ 150; %s The PBS dose rate is 150 Gy/s
Period     = PulseWidth .* 1.1 ; %Period (s)
Meas       =         [5  , 28 , 42 , 50 , 42 ] .* 100 ./ 54; %Skin toxicity score 2.5 Figure 4A
SD         =         [0  , 0  , 0  , 0  , 0  ] ; %Error bars
%Skin toxicity for mice: score 2.5
%Sørensen, B. S., Kanouta, E., Ankjærgaard, C. & Kristensen, L. Proton FLASH : Impact of Dose Rate and Split Dose on Acute Skin Toxicity in a Murine Model. Int. J. Radiat. Oncol. Biol. Phys. 000, 1–11 (2024).

%Compute the [LOOH]f for all the experimental conditions
%--------------------------------------------------------
[LOOHf , AvDoseRate , PkDoseRate ] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2);

LOOHf = NbDelivery .* NbDelivery;
      %We assume that the [LOOH] of each delivery remains stable at the end of the delivery
      %and that the [LOOH] of each delivery adds up to the one from the previous delivery

figure(200)
semilogx(AvDoseRate,LOOHf , '*k');
hold on
xlabel('Average dose rate (Gy/s)')
ylabel('[LOOH]_f (\mu M)')
title(['[LOOH]_f' ])
grid on
drawnow

figure(201)
errorbar(LOOHf , Meas , SD , '*k')
hold on
xlabel('[LOOH]_f (\mu M)')
ylabel('Skin toxicity Score 2.5 (%)')
%title(['AUC' ])
grid on
hold on
for idx = 1:numel(AvDoseRate)
  text(LOOHf(idx)+0.5 , Meas(idx), num2str(NbDelivery(idx)))
end


x = LOOHf(1:end-1)'
x = [x , ones(size(x,1),1)]
y = Meas(1:end-1)';
b = x\y %Make a linear regression of [LOOH]f vs log(DR)
LOOHm = b'*x';
Rsq = -abs(1 - sum((y - LOOHm').^2) / sum((y - mean(y)).^2));
fprintf('Coefficient of determination, : %f \n',Rsq)

figure(201)
hold on
plot(x(:,1) , LOOHm , '-r')
