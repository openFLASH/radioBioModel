%% oxygenKinetics
% Compute the concentration of different radical species as a function of time
% for a continuous beam
%
%% Syntax
% |oxygenKinetics|
%
%
%% Description
% |oxygenKinetics| Start computation and display graphs
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
%[9] Chen, X., Zhong, Z., Xu, Z., Chen, L., & Wang, Y. (2010). 2′,7′-Dichlorodihydrofluorescein as a fluorescent probe for reactive oxygen species measurement: Forty years of application and controversy. Free Radical Research, 44(6), 587–604. https://doi.org/10.3109/10715761003709802
%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

clear
close all

store_y = {};
store_INTRooR = {};
legendSTR = [];
legendAUC = [];

%Define simulation parameters
TimeSimuMax = 150; %s %Duration of the simulation
TimeStep = 1e-3; %Step between time points in simulation
doseRate = [1e7 , 1e6, 1e5 , 1e4, 1e3, 500, 100, 75, 50, 25 , 10 ,  1, 0.5 , 0.25 , 0.03]; %Gy/s %average dose rate
TotalDose = [1, 2, 5 , 7, 10 , 15 , 20 , 30];% Gy
O2 = [10,25, 50, 75, 100,200]; % [u-mol/L]. This enables to sweep calculations through different concentrations


param = paramRadiolytic();
fprintf('conversion factor kr-ge = %g \n',ge2kr(1));
param.pH = 7; %pH of the extra vascular tissue. NB: assume buffered solution [H+] and [OH-] are constant
param.R = @constantBeam;
TissueParam.O2 = O2;

%on/off biological reactions
biologyFLAG = true;


ROOrInt=[];
BeamTime = [];
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
          opts = odeset('NonNegative',1:length(y0),'InitialStep',1e-7,'MaxStep',param.t_on,'Jacobian',@JacRadiolysisKinetics);
        else
          fprintf('Not a pulsed beam -- Default ode time step \n')
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

        markers = {'-','--',':','-.'};
        colour = {'k','b','g','r','c','m','y'};

        for fig = 1:length(y0)
          figure(fig)
          %semilogx(t,y(:,fig), markers{doseRateIndex},'MarkerIndices',1:30:length(t),'MarkerSize',5)
          loglog(t,y(:,fig), markers{mod(doseRateIndex,length(markers))+1},'Linewidth', 1.3,'MarkerIndices',1:45:length(t),'MarkerSize',7)
          set(gca,'FontSize',20)
          hold on

          xlabel('Time (s)', 'FontSize', 24)
          ylabel('Concentration (\mu mol/l)', 'FontSize', 24)
          title(['[',labels{fig},'] evolution',], 'FontSize', 24)
          legendSTR{fig,legendIndex} = ['dD/dt = ',num2str(param.R0,'%1.2g'),'Gy/s'];
          leg = legend(legendSTR{fig,:});
          set(leg,'FontSize',14);

          if(setAxis)
            axis([TimeStep/10,tspan(end),1e-3.*max(max(y(:,fig))),1.1.*max(max(y(:,fig)))])
            setAxis = 0;
          end
          grid on
        end

        legendIndex = legendIndex+1;

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


        % Store INTRooR per O2/totalDose
        store_INTRooR{iO2} = INTRooR;
        store_RooR{doseRateIndex, doseIndex,iO2} = RooR;


        %Compute integral under H2O2(t) curve
        %=====================================
        H2O2 = y(:,3);
        INTH2O2(doseRateIndex, doseIndex) = trapz(t,H2O2).*1e3; %nM.s
        fprintf('Integral H2O2 = %g uM.s \n',INTH2O2(doseRateIndex, doseIndex))

        drawnow
        fprintf('DONE \n')
    end %for doseRateIndex
  end % for doseIndex
end %for iO2

AUC_file = 'AUC_ROO_Dstudy_k.mat';
doseRate_file = 'doseRate_Dstudy_k.mat';
TotalDose_file = 'TotalDose_Dstudy_k.mat';
O2_file = 'O2_file_Dstudy_k.mat';
ROOr_file = 'ROOr_Dstudy_k.mat'

save(AUC_file,'store_INTRooR')
save(doseRate_file,'doseRate')
save(TotalDose_file,'TotalDose')
save(O2_file,'O2')
save(ROOr_file ,'store_RooR')


data = store_y{1}{1}{1};
OHr = data(:,5) .* 1e-6; %u-mol/l
t = data(:,1);
figure(100)
loglog(t,OHr)


%==================================
% system of Ordinary differential Equation
% Kinetics of oxygen depletion by aqueous electrons
%==================================

function dydt = odefcnComprehensive(t,y,param,TissueParam,biologyFLAG)
  dydt = radiolysisKinetics(t,y,param,TissueParam,biologyFLAG);
end
