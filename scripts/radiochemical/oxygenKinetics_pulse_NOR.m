% Compute the concentration [LOOH]f
% Use the radio-kinetic model for the 2 phases model
% for the experimental conditions of the Novel Object Recognition test
%
%% Syntax
%
%REFERENCE
% Montay-Gruel, P. et al. Irradiation in a flash: Unique sparing of memory in mice after whole brain irradiation with dose rates above 100 Gy/s. Radiother. Oncol. 124, 365–369 (2017).
% Montay-gruel, P., Acharya, M. M., Petersson, K., Alikhani, L. & Yakkala, C. Long-term neurocognitive benefits of FLASH radiotherapy driven by reduced reactive oxygen species. PNAS Latest Artic. (2019). doi:10.1073/pnas.1901777116
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

TotalDose  = [10       , 10  , 10  , 10   , 10    , 10    , 10    , 10    , 10    , 10    , 10  ];% Gy
O2         = [100      , 50  , 50  , 50   , 50    , 50    , 50    , 50    , 50    , 50    , 50  ]; %u-mol/l
NbPulses   = [1        , 1000, 1000, 333  , 100   , 50    , 33    , 17    , 10    , 2     ,  1  ];
PulseWidth = [2        , 2   , 2   , 2    , 2     , 2     , 2     , 2     , 2     , 2     , 2   ] .* 1e-6 ; %s
Period     = [2.1e-3   , 100 , 10  , 10   , 10    , 10    , 10    , 10    , 10    , 10    , 2.1e-3] .* 1e-3 ; %Period (s)
Meas       = [64       , 53  ,  56 , 56.7 , 57.1  , 54.4  , 65.9  , 72.1  , 75.9  , 79.1  , 75.9]; %Recognition ratio, figure 1
SD         = [2        , 0.8 , 1.2 , 2.7  , 1.5   , 1.9   ,  4    , 1.4   , 1.8   , 1.4   , 2   ] .* (50./11.6) ./ 2; %Error bars
doseRateEXP= [5.60E+06 , 0.1 ,  1  , 3    , 10    , 20    , 30    , 60    , 100   , 500   , 5.60E+06]; %Gy/s
% Novel Object Recognition test
% Data point 2 to end:
% Montay-Gruel, P. et al. Irradiation in a flash: Unique sparing of memory in mice after whole brain irradiation with dose rates above 100 Gy/s. Radiother. Oncol. 124, 365–369 (2017).
% First data points : with carbogen
% Montay-gruel, P., Acharya, M. M., Petersson, K., Alikhani, L. & Yakkala, C. Long-term neurocognitive benefits of FLASH radiotherapy driven by reduced reactive oxygen species. PNAS Latest Artic. (2019). doi:10.1073/pnas.1901777116


%Compute the [LOOH]f for all the experimental conditions
%--------------------------------------------------------
[LOOHf , AvDoseRate , PkDoseRate ] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2);


figure(200)
semilogx(AvDoseRate(2:end),LOOHf(2:end) , '*k','MarkerSize',15);
hold on
semilogx(AvDoseRate(1),LOOHf(1) , '*r','MarkerSize',15);
xlabel('Average dose rate (Gy/s)', 'FontSize', 24)
ylabel('[LOOH]_f (\mu M)', 'FontSize', 24)
title(['[LOOH]_f' ], 'FontSize', 24)
grid on
drawnow

figure(201)
errorbar(LOOHf(2:end-1) , Meas(2:end-1) , SD(2:end-1) , '*k','MarkerSize',15)
set(gca,'FontSize',20)
hold on
errorbar(LOOHf(1) , Meas(1) , SD(1) , 'or','MarkerSize',15)
errorbar(LOOHf(end) , Meas(end) , SD(end) , '^g','MarkerSize',15)
xlabel('[LOOH]_f (\mu M)', 'FontSize', 24)
ylabel('Recognition ratio (%)', 'FontSize', 24)

grid on
% hold on
% for idx = 1:numel(AvDoseRate)
%   text(LOOHf(idx)+0.1 , Meas(idx), ['DRa=' num2str(AvDoseRate(idx),'%3.2g'), ' Gy/s - DRp=' num2str(PkDoseRate(idx),'%3.2g'), ' Gy/s'])
% end


x = LOOHf(1:end-2)'
x = [x , ones(size(x,1),1)]
y = Meas(1:end-2)';
b = x\y %Make a linear regression of [LOOH]f vs log(DR)
LOOHm = b'*x';
Rsq = -abs(1 - sum((y - LOOHm').^2) / sum((y - mean(y)).^2));
fprintf('Coefficient of determination, : %f \n',Rsq)

figure(201)
hold on
plot(x(:,1) , LOOHm , '-r','Linewidth', 1.3)
