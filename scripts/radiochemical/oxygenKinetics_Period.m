%% oxygenKinetics_Period
% Compute the concentration of different radical species as a function of time
% for a pulsed beam
%
%% Syntax
% |oxygenKinetics_Period|
%
%
%% Description
% |oxygenKinetics_Period| Start computation and display graphs
%
%
%% Input arguments
% None
%
%
%% Output arguments
%
% None
%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

clear
close all

legendSTR = [];
legendAUC = [];

%Define simulation parameters
TimeSimuMax = 150; %s %Duration of the simulation
TimeStep = 1e-3; %Step between time points in simulation
doseRate = 50; %Gy/s
TotalDose = 10; %Gy
O2 = 50; %u-mol/l
Period    = [100, 50, 20  ].*1e-3;%s
DutyCycle = [0.5, 1e-1 , 1e-2]; %fraction

[param , TissueParam ] = paramKinectisDiffusion();
fprintf('conversion factor kr-ge = %g \n',ge2kr(1));
param.pH = 7; %pH of the extra vascular tissue. NB: assume buffered solution [H+] and [OH-] are constant

BeamFunction = {'pulsedBeam'};

beamType = 1;
doseRateIndex =1;
doseIndex =1;
iO2 = 1;
ROOrInt=[];
BeamTime = [];

%on/off biological reactions
biologyFLAG = true;
if (biologyFLAG)
    fprintf('Biological reactions on \n');
else
    fprintf('Biological reactions off \n');
end

indexLegend = 1;
INTRooR = [];
INTH2O2 = [];
INTRr = [];
legendIndex = 1;

for PeriodIndex = 1:length(Period)
for DutyCycleIndex = 1:length(DutyCycle)

    param.R = str2func(BeamFunction{beamType});
    param.R0 = doseRate(doseRateIndex); %Average dose rate
    param.T = Period(PeriodIndex); %s
    param.t_on = param.T .* DutyCycle(DutyCycleIndex); %s


    fprintf('Computing for average dose rate %g Gy/s .... \n',param.R0);
    fprintf('Period  : %g s \n',param.T)
    fprintf('Beam ON : %g s \n',param.t_on)
    Rav = param.R0;

    % use if continuous irradiation
    param.td = TotalDose(doseIndex) ./ Rav; % s Beam ON time. Set the beam time to deliver the same total dose for all dose rate
    BeamTime(doseRateIndex) = param.td;
    fprintf('Computing for dose %g Gy => delivery time %g s \n',TotalDose(doseIndex),BeamTime(doseRateIndex));

    TimeSimuMaxLocal = BeamTime(doseRateIndex).*1.1; %the simulation is run only during the beam ON time
    tspan =  0:param.t_on./10:TimeSimuMaxLocal;
    tspanmax = max(tspan);
    TimeSimuMaxLocal = TimeSimuMax;
    tspan = [tspan , (tspanmax + param.t_on./10):1:TimeSimuMaxLocal];
    fprintf('Simulation end time = %f s \n',tspan(end))


    %Initial concentrations
    %=======================
    [~ , labels]= radiolysisKinetics();
    Nbconcentrations = length(labels);
    fprintf('Number of tracked concentrations: %d \n',Nbconcentrations);
    labels

    y0 = zeros(1,Nbconcentrations); %Set all the radicals concentrations to zero
    y0(2)= O2(iO2); %u-mol/l Oxygen inital concentration
    fprintf('[O_2]_{ini} : %g uM \n',O2(iO2))

    if (strcmp(func2str(param.R) , 'pulsedBeam'))
      fprintf('Pulsed beam -- Max step = %g s \n',param.t_on);
      opts = odeset('NonNegative',1:length(y0),'InitialStep',1e-8,'MaxStep',param.t_on,'OutputFcn','odeprint');
    else
      fprintf('Not a pulsed beam -- Default ode time step \n') %Refined max step 1e-5
      opts = odeset('NonNegative',1:length(y0));
    end
    Tstart = datetime;
    fprintf('Computation starts at %s \n',Tstart);
    [t,y] = ode15s(@(t,y) odefcnComprehensive(t,y,param, TissueParam, biologyFLAG),tspan,y0,opts);
    Tend = datetime;
    fprintf('Computation ends at %s \n',Tend);
    fprintf('Duration : %s \n',Tend-Tstart);

    [~ , Rp] = param.R(0,param);
    fprintf('Peak dose rate = %f Gy/s \n',Rp);

    %Display graph for all species
    %=============================

    setAxis = 1;
    [~ , labels]= radiolysisKinetics();
    markers = {'-+','-o','-v','-s','-p','-h','-d','-^','-*','-+','-o','-v','-s','-p','-h','-d','-^','-*'};
    colour={'k','b','r','g','c','m','y'};
    for fig = 1:length(y0)
      figure(fig)
      semilogx(t,y(:,fig), [markers{mod(DutyCycleIndex,length(markers))+1},colour{mod(PeriodIndex,length(colour))+1}],'MarkerIndices',1:30:length(t),'MarkerSize',5)
      set(gca,'FontSize',10)
      hold on

      xlabel('Time (s)', 'FontSize', 16)
      ylabel('Concentration (\mu mol/l)', 'FontSize', 16)
      title(['[',labels{fig},'] evolution',], 'FontSize', 16)
      legendSTR{fig,legendIndex} = ['T = ',num2str(param.T,'%1.2g'),'s -- t_{on} : ',num2str(param.t_on,'%1.2g s')];
      legend(legendSTR{fig,:});

      if(setAxis)
        axis([TimeStep/10,tspan(end),1e-3.*max(max(y(:,fig))),1.1.*max(max(y(:,fig)))])
        setAxis = 0;
      end
      grid on
    end


    %Compute integral under Rr(t) curve
    %=====================================
    [INTRrTMP,INTRooRTMP,INTH2O2TMP] = computeIntegrals(t,y);

    INTRr(PeriodIndex, DutyCycleIndex) = INTRrTMP;
    fprintf('Period  : %g s \n',param.T)
    fprintf('Beam ON : %g s \n',param.t_on)
    fprintf('Integral R. = %g uM.s \n',INTRr(PeriodIndex, DutyCycleIndex))

    INTRooR(PeriodIndex, DutyCycleIndex) = INTRooRTMP;
    fprintf('Integral ROO. = %g uM.s \n',INTRooR(PeriodIndex, DutyCycleIndex))

    INTH2O2(PeriodIndex, DutyCycleIndex) = INTH2O2TMP; %nM.s
    fprintf('Integral H2O2 = %g uM.s \n',INTH2O2(PeriodIndex, DutyCycleIndex))

    drawnow
    fprintf('DONE \n')

    legendIndex = legendIndex+1;
  end %for DutyCycleIndex

  figure(1)
  semilogx(DutyCycle.*100 , squeeze(INTRooR(PeriodIndex, :)),'-o');
  xlabel('Duty cycle (%)')
  ylabel('Exposure to ROO^. (uM.s)')
  legendSTR{legendIndex} = ['Freq = ',num2str(1e-3./param.T,'%1.2g'),'kHz'];
  legend(legendSTR{:});
  grid on
  drawnow
  hold on
  %legendIndex = legendIndex+1;

end  %for PeriodIndex


INTRooR



%==================================
% system of Ordinary differential Equation
% Kinetics of oxygen depletion by aqueous electrons
%==================================

function dydt = odefcnComprehensive(t,y,param, TissueParam, biologyFLAG)
    dydt = radiolysisKinetics(t,y,param, TissueParam,biologyFLAG);
end


%======================================================================
function [INTRr , INTRooR , INTH2O2] = computeIntegrals(t,y)

  %Compute integral under Rr(t) curve
  %=====================================
  Rr = y(:,8);
  INTRr = trapz(t,Rr).*1e3; %nM.s


  %Compute integral under ROOr(t) curve
  %=====================================
  RooR = y(:,9);
  INTRooR = trapz(t,RooR).*1e3; %nM.s


  %Compute integral under H2O2(t) curve
  %=====================================
  H2O2 = y(:,3);
  INTH2O2 = trapz(t,H2O2).*1e3; %nM.s

end
