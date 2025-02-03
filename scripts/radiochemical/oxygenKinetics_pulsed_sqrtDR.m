%---------------------------------------------
% Check that the model predict 1/sqrt(DR) dependence of [LOOHf]
%	Mihaljevich, Branka and Tartaro, Ivana and Ferreri, Carla and Chatgilialoglu, Chryssostomos.  Linoleic acid peroxidation vs. isomerization: a biomimetic model of free radical reactivity in the presence of thiols.  Org. Biomol. Chem., 9:3541-3548, 2011.
%---------------------------------------------
clear
close all

colors = {'k','b','g','r','c','m','y'};

colors = {'b','r'};
symbols = {'_', '^' , 'o' };

TotalDose =  10  ;% Gy
Period = [1e-3 , 10e-3]; %Period
DutyCycle = [0.1 , 0.5 , 0.9];
NbPulses = [1 , 2 , 10 , 12 , 15 , 25 , 50 ,1e2 , 1e3 ]; %Make sure that there is an integral number of complete pulses
O2 = 50; %u-mol/l


for Ti = 1:numel(Period)
  for Qi = 1:numel(DutyCycle)

    DRa(Ti,Qi,:) =  TotalDose ./ (NbPulses .* Period(Ti)); %Gy/s %average dose rate
    PulseWidth = DutyCycle(Qi) .*  Period(Ti);

    for Ni = 1:numel(NbPulses)
        LOOHf(Ti,Qi,Ni)  = getLOOHf(TotalDose , Period(Ti) , PulseWidth , NbPulses(Ni) , O2 , [] , false);
    end

    %Plot  1 / sqrt(DR_a)
    figure(11)
    plot(squeeze(1./sqrt(DRa(Ti,Qi,:))) , squeeze(LOOHf(Ti,Qi,:)) , [symbols{Qi} colors{Ti}])
    hold on
    drawnow

  end
end

figure(11)
xlabel('$1/\sqrt{DR_a} (Gy/s)^{-2}$ ',Interpreter='latex')
ylabel('[LOOH]_f \muM')
grid


%Fit a line on the sqrt(DR) plot < 0.2 (Gy/s)-2
%---------------------------------------
xt = 1./sqrt(DRa(:));
yt = LOOHf(:);
wF = xt > 0.2; %1/sqrt(DRa) > 0.2
xt(wF) = [];
yt(wF) = [];

xt = [xt , ones(size(xt,1),1)];

b = xt\yt %Make a linear regression of [LOOH]f vs log(DR)

LOOHm = b'*xt';
Rsq = 1 - sum((yt - LOOHm').^2) / sum((yt - mean(yt)).^2);
fprintf('Coefficient of determination, : %f \n',Rsq)

%plot ther fitted line on the graph
x = 0:0.05:0.3;
x = [x' , ones(size(x',1),1)];
LOOHm = b'*x';

figure(11)
hold on
plot(x(:,1),LOOHm,'-g')






function legendSTR = addcurve(Name , t,y,labels, legendSTR )
  idx = find(strcmp(labels , Name));
  loglog(t,y(:,idx))
  legendSTR{end+1} = Name;
end
