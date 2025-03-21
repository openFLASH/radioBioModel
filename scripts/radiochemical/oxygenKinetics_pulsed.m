%------------------------------------
% PLot the figures of the article on the pulsed beam with 2 phases
%------------------------------------

clear
close all

colors = {'k','b','g','r','c','m','y'};
symbols = {'o', '*', '+',  '.', 'x', '_', '|', 'square'	, 'diamond'	, '^'	, 'v'	, '>'	, '<'	, 'pentagram', 'hexagram'};


% Weiss
% Weiss H, Epp ER, Heslin JM, Ling CC, Santomasso A. Oxygen depletion in cells irradiated at ultra-high dose-rates and at conventional dose-rates. Int J Radiat Biol 1974;26:17–29
%------------------------
TotalDose  = 10; %Gy
Period  = 0.2; %s
PulseWidth = 0.2; %s
NbPulses  = 1
O2 = 50; %u-mol/l

[~ , ~ , ~ , ~ , t , y,labels] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , [] , true);
%[t, y , labels] = getConc(TotalDose , Period , PulseWidth , NbPulses , O2 , [] , true);


figure(100)
legendSTR = {};
legendSTR = addcurve('LOO^.' , t,y,labels, legendSTR );
hold on
legendSTR = addcurve('LOOH' , t,y,labels, legendSTR );
legendSTR = addcurve('e^-_{aq}' , t,y,labels, legendSTR );

xlabel('Time (s)')
ylabel('Concentration (\mu mol/l)')
legend(legendSTR)
grid minor


figure(101)
legendSTR = {};
legendSTR = addcurve('O_2' , t,y,labels, legendSTR );

xlabel('Time (s)')
ylabel('Concentration (\mu mol/l)')
title('[O_2] evolution')
grid minor



%Howard Flanders
%	Howard-Flanders P, Moore D. The time interval after pulsed irradiation within which injury to bacteria can be modified by dissolved oxygen. I. A search for an effect of oxygen 0.02 second after pulsed irradiation. Radiat Res 1958;9:422–37.
%------------------------
TotalDose  = 100; %Gy
Period  = 90; %s
PulseWidth = 90; %s
NbPulses  = 1
O2 = 0; %u-mol/l

%[t, y , labels] = getConc (TotalDose , Period , PulseWidth , NbPulses , O2 , [] , true);
[~ , ~ , ~ , ~ , t , y,labels] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , [] , true);

figure(102)
legendSTR = {};
legendSTR = addcurve('L^.' , t,y,labels, legendSTR );
hold on
legendSTR = addcurve('e^-_{aq}' , t,y,labels, legendSTR );

xlabel('Time (s)')
ylabel('Concentration (\mu mol/l)')
title('[L^.] evolution')
grid minor

%Compare AUC and [LOOH]f
%------------------------

TotalDose  = 10; %Gy
Period  = 1e-3; %s
PulseWidth = [2,5,10,50,100,500,900] .* 1e-6; %s
NbPulses  = [2 , 1e1 , 5e1 , 1e2 , 5e2 , 1e3 , 5e3 ];
O2 = 50; %u-mol/l

legendSTR = {};

for pw = 1:numel(PulseWidth)
    for nbp = 1:numel(NbPulses)
        %[t, y , labels] = getConc(TotalDose , Period , PulseWidth(pw) , NbPulses(nbp) , O2 , [] , false);
        [~ , ~ , ~ , ~ , t , y,labels] = getLOOHf(TotalDose , Period , PulseWidth(pw) , NbPulses(nbp) , O2 , [] , false);
        idx = find(strcmp(labels , 'LOOH'));
        LOOHf(nbp , pw) =  y(end,idx); %uM
        idx = find(strcmp(labels , 'LOO^.'));
        INTRooR(nbp , pw) = trapz(t,y(:,idx)).*1e3; %nM.s
        DRa(nbp) =  TotalDose ./ (NbPulses(nbp) .* Period); %Gy/s %average dose rate
    end

    figure(110)
    semilogx(DRa , squeeze(LOOHf(: , pw)) , ['-' symbols{mod(pw,length(symbols))+1} colors{mod(pw,length(colors))+1}])
    hold on
    title('[LOOH]_f')
    legendSTR{end+1} = ['Pulse width = ' num2str(PulseWidth(pw).*1e6) ' \mus'];
    drawnow

    figure(111)
    semilogx(DRa , squeeze(INTRooR(: , pw)) , ['-' symbols{mod(pw,length(symbols))+1} colors{mod(pw,length(colors))+1}])
    hold on
    title('AUC')
    drawnow

end

figure(110)
grid minor
legend(legendSTR)
xlabel('Average dose rate (Gy/s)')
ylabel('[LOOH]_f (\mu M)')


figure(111)
grid minor
legend(legendSTR)
xlabel('Average dose rate (Gy/s)')
ylabel('AUC (nM.s)')

figure(112)
for pw = 1:numel(PulseWidth)
  semilogx(DRa , 1e-3.*squeeze(INTRooR(: , pw))./squeeze(LOOHf(: , pw)) , ['-' symbols{mod(pw,length(symbols))+1} colors{mod(pw,length(colors))+1}])
  hold on
end
grid minor
legend(legendSTR)
xlabel('Average dose rate (Gy/s)')
ylabel('AUC / [LOOH]_f (s)')

%------------------------
% Add the curve of concentration vs time to a plot
%------------------------
function legendSTR = addcurve(Name , t,y,labels, legendSTR )
  idx = find(strcmp(labels , Name));
  loglog(t,y(:,idx))
  legendSTR{end+1} = Name;
end
