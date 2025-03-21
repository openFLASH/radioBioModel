% Compute the concentration [LOOH]f
% Use the radio-kinetic model for the 2 phases model
% for the experimental conditions of cyclotron pulse structures
%
%% Syntax
%
%REFERENCE
%Reference: Karsch, L. et al. Beam pulse structure and dose rate as determinants for the flash effect observed in zebrafish embryo. Radiother. Oncol. 173, 49–54 (2022).
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

TotalDose =  31  ;% Gy

O2         = [17     ,  17         , 17             ,  17         , 17       , 17]; %u-mol/l
NbPulses   = [1      ,  1          , 5              ,  1          ,  1       , 1]; %Make sure that there is an integral number of complete pulses
PulseWidth = [239    ,  0.09       , 60e-6          , 0.29e-3     , 199      , 0.1] ; %s
Period     = [240    ,  0.1        , 40e-3          , 0.3e-3      , 200      , 0.115] ; %Period (s)
Label      = {'CONV' ,'UHDR_{iso}' ,'UHDR_{synchro}','UHDR_{max}' , 'P_{ref}','P_{UHDR}'};
Meas       = [3402.7 , 3525.6      , 3517.4         , 3636.7      , 3505.1   , 3643.3]; %um Mean body length of zebra fish. Table 2
SD         = [120.8  , 98.5        , 110.2          , 106.5       , 150.4    , 84.9 ]; %Error bars
%Reference: Karsch, L. et al. Beam pulse structure and dose rate as determinants for the flash effect observed in zebrafish embryo. Radiother. Oncol. 173, 49–54 (2022).
% 'CONV' ,'UHDR_{iso}' and 'UHDR_{max}' have the 13 MHz modulation => consider quasi continuous => DRp = DRa

%Save the simulation parameters
data.TotalDose = TotalDose;
data.Period =Period ;
data.PulseWidth = PulseWidth ;
data.NbPulses = NbPulses;
data.O2 = O2;
save(fullfile(folder , ['data_' num2str(Period) '.mat']) , 'data')


%Rate constants
%------------------
% [kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants();
%
% %Manually change the value of some rate constants in the model
% param = getDefaultParam();
% kName = {  'kbr2' , 'kb3' , 'kLOOself' , 'kb8' , 'kbr'};
% kValue = [kbr2 , kb3,  kLOOself , kb8 , kbr];
% param = set_k_in_param(param,kValue,kName); %Update the value of the rate constant
%
% AvDoseRate =  TotalDose ./ (NbPulses .* Period); %Gy/s %average dose rate

%Loop for every simulation
for idx = 1:numel(O2)


    %on/off biological reactions
    indexLegend = 1;
    legendSTR = {};
    legendSTR2 = {};


    [LOOHf(idx) , AvDoseRate(idx) , ~ , ~ , t , y] = getLOOHf(TotalDose , Period(idx) , PulseWidth(idx) , NbPulses(idx) ,  O2(idx) );

    % param.R0   = AvDoseRate(idx); %Average dose rate
    % param.T    = Period(idx); %s
    % param.t_on = PulseWidth(idx); %s Duration of a single pulse
    %
    % legendSTR{end+1} =  ['<dD/dt> = ' num2str(param.R0,'%2.1g') ' Gy/s -- T_{pulse} = ' num2str(param.t_on .* 1e3,'%2.1g') ' ms'];
    % fprintf('Computing for average dose rate %1.3g Gy/s .... \n',param.R0);
    % fprintf('Period  : %g s \n',param.T)
    % fprintf('Beam ON : %1.3g s \n',param.t_on)
    %
    % param.td = TotalDose ./ param.R0; % s Beam ON time. Set the beam time to deliver the same total dose for all dose rate
    % fprintf('Computing for dose %g Gy => delivery time %g s \n',TotalDose,param.td);
    % fprintf('Number of pulses : %d \n',ceil(param.td ./ param.T))
    %
    % %Initial concentrations
    % %=======================
    % [dydt , labels]= radiolysisKinetics2P_a();
    % Nbconcentrations = length(labels);
    % fprintf('Number of tracked concentrations: %d \n',Nbconcentrations);
    % labels
    %
    % [y0 , t0] = getY0(param , @radiolysisKinetics2P_a); %Radical concentration at begining of homogeneous chemical phase
    % y0(2)  = O2(idx); %u-mol/l Oxygen inital concentration
    %
    % if (strcmp(func2str(param.R) , 'pulsedBeam'))
    %   fprintf('Pulsed beam -- Pulse width = %g s \n',param.t_on);
    %   opts = odeset('NonNegative',1:length(y0),'InitialStep',param.t_on.*1e-4,'Jacobian',@JacRadiolysisKinetics2P_a);
    % else
    %   fprintf('Not a pulsed beam -- Default ode time step \n')
    %   opts = odeset('NonNegative',1:length(y0),'Jacobian',@JacRadiolysisKinetics2P_a);
    % end
    % Tstart = datetime;
    % fprintf('Computation starts at %s \n',Tstart);
    % [t,y] = pulsedODE(@odefcnComprehensive, t0 , y0,opts,param,param);
    % Tend = datetime;
    % fprintf('Computation ends at %s \n',Tend);
    % fprintf('Duration : %s \n',Tend-Tstart);

    % Store all computation results per O2/totalDose/doseRate
    store_y{idx} = [t,y];

    % [~ , Rp] = param.R(0,param);
    % fprintf('Peak dose rate = %2.3g Gy/s \n',Rp);
    %
    %Display graph for all species
    %=============================
    displayGraphSpecies(t, y , idx , legendSTR , idx);

    %Get final LOOH concentration
    %=====================================
    %LOOHf(idx) =  y(end,10); %uM
    % fprintf('Average dose rate = %f Gy/s \n',param.R0)
    % fprintf('Dose %g Gy => delivery time %g s \n',TotalDose,param.td);
    % fprintf('[LOOH]f = %g uM \n',LOOHf(idx))

    %Compute integral under ROOr(t) curve
    %=====================================
    RooR = y(:,9);
    INTRooR(idx) = trapz(t,RooR).*1e3; %nM.s

    fprintf('DONE \n')

end %for idx


figure(200)
semilogx(AvDoseRate,LOOHf , '*k');
hold on
xlabel('Average dose rate (Gy/s)')
ylabel('[LOOH]_f (\mu M)')
title(['[LOOH]_f' ])
grid on
drawnow


figure(201)
semilogx(AvDoseRate,INTRooR, '*k');
hold on
xlabel('Average dose rate (Gy/s)')
ylabel('AUC (nM.s)')
title(['AUC' ])
grid on
drawnow

figure(202)
semilogx(AvDoseRate,1e-3 .* INTRooR ./ LOOHf, '*k');
hold on
xlabel('Average dose rate (Gy/s)')
ylabel('AUC ./ [LOOH]_f (s) ')
title(['AUC' ])
grid on
drawnow

figure(203)
errorbar(LOOHf(1:4) , Meas(1:4) , SD(1:4), '*k')
xlabel('[LOOH]_f (\mu M)')
ylabel('Zebra fish length (\mum) ')
%title(['AUC' ])
grid on
drawnow


x = LOOHf(1:4)'
x = [x , ones(size(x,1),1)]
y = Meas(1:4)';
b = x\y %Make a linear regression of [LOOH]f vs log(DR)
LOOHm = b'*x';
Rsq = -abs(1 - sum((y - LOOHm').^2) / sum((y - mean(y)).^2));
fprintf('Coefficient of determination, : %f \n',Rsq)

figure(203)
hold on
plot(LOOHf(1:4) , LOOHm , '-r')
for idx = 1:4
  text(LOOHf(idx) , Meas(idx),Label{idx})
end
