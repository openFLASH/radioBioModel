%% oxygenKinetics_kSensitivity
% Compute the concentration of different radical species as a function of time
% Vary the value of the rate constants
%
%% Syntax
% |oxygenKinetics_kSensitivity|
%
%
%% Description
% |oxygenKinetics_kSensitivity| Start computation and display graphs
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
%REFERENCE
%[1] Azzam, E. I., Jay-Gerin, J.-P., Pain, D. (2012). Ionizing radiation-induced metabolic oxidative stress and prolonged cell injury. Cancer Letters, 327(1–2), 48–60. http://doi.org/10.1016/j.canlet.2011.12.012
%[2] Buxton, G. V, Greenstock, C. L., Helman, W. P., Ross, A. B. (1988). Critical Review of Rate Constants for Reactions of Hydrated Electrons, Hydrogen Atoms and Hydroxyl Radicals (OH/O−) in Aqueous Solution. Atomic Energy, 17, 513–886. http://doi.org/10.1063/1.555805
%[3] Free Radicals in Biology, Volume 3, Ed Willimam A. Pryor, Academic press 1977.
%[4] Cobut, V., Frongillo, Y., Patau, J. P., Goulet, T., Fraser, M. J., Jay-Gerin, J. P. (1998). Monte Carlo simulation of fast electron and proton tracks in liquid water - I. Physical and physicochemical aspects. Radiation Physics and Chemistry, 51(3), 229–243
%[5] Frongillo, Goulet, Fraser, Cobut, Patau and Jay-Gerin Monte carlo simulation of fast electron and proton tracks in liquid water ii. nonhomogeneous chemistry radiat. phys. chem. vol. 51, no. 3, pp. 245-254, 1998
%[6] Kreipl et al. Time- and space-resolved Monte Carlo study of water radiolysis for photon, electron and ion irradiation Radiat Environ Biophys (2009) 48:11–20
%[7] Introduction: radiolysis. (n.d.). http://doi.org/10.1007/s00024-009-0507-0
%[8] Weiss, H. (1972). An Equation for Predicting the Surviving Fraction of Cells Irradiated with Single Pulses Delivered at Ultra-High Dose Rates. Radiation Research, 50(2), 441–452. Retrieved from https://www.jstor.org/stable/3573501?seq=1#metadata_info_tab_contents%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

clear
close all

store_y = {};
store_INTRooR = {};

legendSTR = [];
legendAUC = [];
ROOrInt=[];
BeamTime = [];

TimeSimuMax = 150; %s %Duration of the simulation
TimeStep = 1e-3; %Step between time points in simulation
simu = 'dose'; %The dose is varied, [O2] is constant
%simu = 'O2'; %The dose is constant, [O2] is varied

switch simu
case 'dose'
  doseRate = [1e7 , 1e6, 1e5 , 1e4, 1e3, 500, 100, 75, 50, 25 , 10 , 9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5 , 0.25 , 0.03]; %Gy/s %average dose rate
  TotalDose = [1, 2, 5 , 7, 10 , 15 , 20 , 30 ];% Gy
  O2 = 50; %u-mol/l
case 'O2'
  doseRate = [1e7 , 1e6, 1e5 , 1e4, 1e3, 500, 100, 75, 50, 25 , 10 , 9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5 , 0.25 , 0.03]; %Gy/s %average dose rate
  O2 = [10,25, 50, 75, 100,200]; % [u-mol/L]. This enables to sweep calculations through different concentrations
  TotalDose = 10; %Gy
end

param = paramRadiolytic();
TissueParam.O2 = O2;

%Manually change the value of some rate constants in the model
kName = { 'kb2'           , 'kbr'    , 'kbr2'           , 'kb3'       , 'kb11'    , 'kb8'    , 'kROOself'     , 'AUCmax'  }
k0    = [  3.2388222481e4 , 3.924299 , 7.202635958249e6 ,  712.793261 , 91.241443 , 0.155211 , 1.2791236753e4 ,  25.457972];
param = set_k_in_param(param,k0,kName);

fprintf('conversion factor kr-ge = %g \n',ge2kr(1));
param.pH = 7; %pH of the extra vascular tissue. NB: assume buffered solution [H+] and [OH-] are constant

%param.R = @pulsedBeam;
param.R = @constantBeam;

%on/off biological reactions
biologyFLAG = true;
if (biologyFLAG)
    fprintf('Biological reactions on \n');
else
    fprintf('Biological reactions off \n');
end

indexLegend = 1;

for iO2 = 1:length(O2)
  INTRooR = [];
  INTH2O2 = [];
  INTRr = [];
  legendIndex =1;

  for doseIndex = 1:length(TotalDose)

    for doseRateIndex = 1:length(doseRate)
        param.R0 = doseRate(doseRateIndex); %Average dose rate
        fprintf('Computing for average dose rate %g Gy/s .... \n',param.R0);
        Rav = param.R0;

        % use if continuous irradiation
        param.td = TotalDose(doseIndex) ./ Rav; % s Beam ON time. Set the beam time to deliver the same total dose for all dose rate
        BeamTime(doseRateIndex) = param.td;
        fprintf('Computing for dose %g Gy => delivery time %g s \n',TotalDose(doseIndex),BeamTime(doseRateIndex));

        if (BeamTime(doseRateIndex) > TimeSimuMax)
            TimeSimuMaxLocal = 1.2.*TotalDose(doseIndex)./doseRate(doseRateIndex); %s Total duration of the simulation
            warning(['Beam ON time is beyond end of simulation. Adapting. Max sim time = ',num2str(TimeSimuMaxLocal,'%1.2g'),' s']);
        else
            TimeSimuMaxLocal = TimeSimuMax;
        end

        tlog = -9:0.1:-3;
        tspan = 10.^tlog;
        tspan = [tspan , 2e-3:1e-3:TimeSimuMaxLocal];

        fprintf('Simulation end time = %f s \n',tspan(end))


        %Initial concentrations
        %=======================
        [dydt , labels]= radiolysisKinetics();
        Nbconcentrations = length(labels);
        fprintf('Number of tracked concentrations: %d \n',Nbconcentrations);
        labels

        y0 = zeros(1,Nbconcentrations); %Set all the radicals concentrations to zero
        y0(2)= O2(iO2); %u-mol/l Oxygen inital concentration

        if (strcmp(func2str(param.R) , 'pulsedBeam'))
          fprintf('Pulsed beam -- Max step = %g s \n',param.t_on);
          %opts = odeset('NonNegative',1:length(y0),'InitialStep',1e-7,'MaxStep',param.t_on,'OutputFcn','odeprint','Jacobian',@JacRadiolysisKinetics);
          opts = odeset('NonNegative',1:length(y0),'InitialStep',1e-7,'MaxStep',param.t_on,'Jacobian',@JacRadiolysisKinetics);
        else
          fprintf('Not a pulsed beam -- Default ode time step \n')
          %opts = odeset('NonNegative',1:length(y0),'OutputFcn','odeprint','Jacobian',@JacRadiolysisKinetics);
          opts = odeset('NonNegative',1:length(y0),'Jacobian',@JacRadiolysisKinetics);
        end
        Tstart = datetime;
        fprintf('Computation starts at %s \n',Tstart);
        [t,y] = ode15s(@odefcnComprehensive,tspan,y0,opts,param,param,biologyFLAG);
        Tend = datetime;
        fprintf('Computation ends at %s \n',Tend);
        fprintf('Duration : %s \n',Tend-Tstart);

        % Store all computation results per O2/totalDose/doseRate
        store_y{iO2}{doseIndex}{doseRateIndex} = [t,y];

        [~ , Rp] = param.R(0,param);
        fprintf('Peak dose rate = %f Gy/s \n',Rp);

        %Display graph for all species
        %=============================
        setAxis = 1;
        [~ , labels]= radiolysisKinetics();

        %Compute integral under Rr(t) curve
        %=====================================
        Rr = y(:,8);
        INTRr(doseRateIndex, doseIndex) = trapz(t,Rr).*1e3; %nM.s
        fprintf('Average dose rate = %f Gy/s \n',param.R0)
        fprintf('Dose %g Gy => delivery time %g s \n',TotalDose(doseIndex),BeamTime(doseRateIndex));
        fprintf('Integral R. = %g uM.s \n',INTRr(doseRateIndex, doseIndex))

        %Compute integral under ROOr(t) curve
        %=====================================
        RooR = y(:,9);
        INTRooR(doseRateIndex, doseIndex) = trapz(t,RooR).*1e3; %nM.s
        fprintf('Integral ROO. = %g uM.s \n',INTRooR(doseRateIndex, doseIndex))

        store_RooR{doseRateIndex, doseIndex,iO2} = RooR;


        %Compute integral under H2O2(t) curve
        %=====================================
        H2O2 = y(:,3);
        INTH2O2(doseRateIndex, doseIndex) = trapz(t,H2O2).*1e3; %nM.s
        fprintf('Integral H2O2 = %g uM.s \n',INTH2O2(doseRateIndex, doseIndex))

        fprintf('DONE \n')
    end %for doseRateIndex

    figure(1)
    semilogx(doseRate,squeeze(INTRooR(:,doseIndex)),'o');
    legendSTR{end+1} = ['D=',num2str(TotalDose(doseIndex)),'Gy [O_2]_i=',num2str(O2(iO2)),'\mu M']
    hold on
    xlabel('Dose rate (Gy/s)')
    ylabel('AUC (\mu M.s)')
    legend(legendSTR)
    grid on
    drawnow

  end % for doseIndex

  %clip negative values
  wneg = find(INTRooR < 0);
  INTRooR(wneg) = 0;

  % Store INTRooR per O2/totalDose
  store_INTRooR{iO2} = INTRooR;
end %for iO2




%save the results
RateConst.k = k0;
RateConst.kName = kName;

switch simu
case 'dose'
  AUC_file = 'AUC_ROO_kstudy.mat';
  doseRate_file = 'doseRate_kstudy.mat';
  TotalDose_file = 'TotalDose_kstudy.mat';
  O2_file = 'O2_file_kstudy.mat';
  ROOr_file = 'ROOr_kstudy.mat'
  k_file = 'k_kstudy.mat'

case 'O2'
  AUC_file = 'AUC_ROO_O2_kstudy.mat';
  doseRate_file = 'doseRate_O2_kstudy.mat';
  TotalDose_file = 'TotalDose_O2_kstudy.mat';
  O2_file = 'O2_file_O2_kstudy.mat';
  ROOr_file = 'ROOr_O2_kstudy.mat'
  k_file = 'k_O2_kstudy.mat'
end

save(AUC_file,'store_INTRooR')
save(doseRate_file,'doseRate')
save(TotalDose_file,'TotalDose')
save(O2_file,'O2')
save(ROOr_file ,'store_RooR')
save(k_file ,'RateConst')


%==================================
% system of Ordinary differential Equation
% Kinetics of oxygen depletion by aqueous electrons
%==================================

function dydt = odefcnComprehensive(t,y,param,TissueParam,biologyFLAG)
  dydt = radiolysisKinetics(t,y,param,TissueParam,biologyFLAG);
end
