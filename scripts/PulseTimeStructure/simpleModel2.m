
clear
close all

fileName = 'C:\Users\rla\IBA Group\Proton Flash Therapy - Documents\Local effect model to describe FLASH\pulse time structure\logK.csv';
folderGraph = 'C:\Users\rla\IBA Group\Proton Flash Therapy - Documents\Local effect model to describe FLASH\pulse time structure\figures';
folderData = 'C:\Users\rla\IBA Group\Proton Flash Therapy - Documents\Local effect model to describe FLASH\pulse time structure\data';

%Experimental data
%-----------------
Excell = 'Max_Moist_desq_percent.xlsx';
T = readtable(fullfile(folderData , Excell))

Label = T.paper;
P_desq = T.x_M_d_;
particle = T.particle; %0= proton; 1= electron
SweepTime = T.deadTimeBtwSpots_pulses_ms_;  %Sweep time between spots (ms)
DwellTime = T.beamOnTimeSpot_pulses_ms_; %Dwell time (ms)
no_spots = T.numberOfSpots; %Number of PBS spots for protons
numberOfPulses = T.numberOfPulses; %numberof pulses for electrons
totalDose = T.Dose_Gy_; %(Gy) dose
btw_field_time = T.time_painting_s_ .* 1e3; %(ms) time between painting
no_repainting = T.nbPainting;
irradiatedArea = T.irradiatedArea_mm2_;
O2 = T.oxygenLevel .* 2 .* 50e-6; % [O2] mol/l
      % in air is label 0.5 => [O2] = 50 uM
      % in oxygen is label 1 => [O2] = 100 uM

LabelList = unique(Label);

%format the dataset
%----------------
for idx = 1:numel(particle)

  switch particle(idx)
    case 0 %PROTON
        %DoseRateCalc(totalDose , no_spots , irradiatedArea , switch_X , beam_on_spot , btw_field_time , no_repainting);
        [DRp{idx} ,t{idx}] = DoseRateCalc(totalDose(idx) , no_spots(idx) , irradiatedArea(idx) , SweepTime(idx) , DwellTime(idx) , btw_field_time(idx) , no_repainting(idx));

    case 1 %Electron
      DRp{idx} =[] ;
      t{idx} = [];

  end
end


%celanup the data
wNaN = find (isnan(SweepTime) + isnan(DwellTime));
for idx = numel(wNaN):-1:1
    Label(wNaN(idx)) = [];
    P_desq(wNaN(idx)) = [];
    particle(wNaN(idx)) = [];
    SweepTime(wNaN(idx)) = [];
    DwellTime(wNaN(idx)) = [];
    no_spots(wNaN(idx)) = [];
    numberOfPulses(wNaN(idx)) = [];
    totalDose(wNaN(idx)) = [];
    btw_field_time(wNaN(idx)) = [];
    no_repainting(wNaN(idx)) = [];
    irradiatedArea(wNaN(idx)) = [];
    O2(wNaN(idx)) = [];
    DRp(wNaN(idx)) = [];
    t(wNaN(idx)) = [];
end


%Inital model parameters
%-----------------------
if 0

  %Read initial from this .m file
  %      gOHr  OH* scv  L* scv   l* crea   L*+L*    L*+O2   LOO* term   gamma    [LOOH]50
  %k0 = [ 1 , 1.0247 , 1.0632 , 0.97211 , 0.9938 , 0.97856 , 0.9686 , 0.25 , 0.01 , 11.97 , 22.9 , 36.585 , 11.287];
  %k0 = [ 1 , 1.2368 , 6.2477 , 0.47087 , 0.8867 , 1.0235 , 1.7871 , 0.44249 , 0.037375 , 0.00015 ];
  %k0 = [ 1 , 1.2368 , 6.2477 , 0.47087 , 0.8867 , 1.0235 , 1.7871 , 0.65e-4 , 6e-6 , 1 ];
  %k0 = [ 1 , 1.2845 , 5.6821 , 0.49031 , 1.0063 , 0.9188 , 1.657 , 6.5e-05 , 6e-06 ,   1 ];
  %k0 = [1 , 1.3487 , 5.6821 , 0.49031 , 1.0063 , 0.9188 , 1.657 , 1 , 10e-6 , 45e-6 ]
  %k0 = [1 , 1.342 , 5.6693 , 0.49155 , 1.039 , 0.91592 , 1.6603 , 1.0202 , 3.304e-06 , 4.0352e-05];

  k0 = [1 , 1.2478 , 4.6569 , 0.55653 , 1.1224 , 0.91571 , 1.8648 , 97.37 , 0.13416 , 1.6575];

else
  %Read initial guess from last line of log file.
  %This restarts the optimisation from where it stopped
  fprintf('Loading k0 from file \n')
  Tk = readtable(fileName);
  NbCol = size(Tk,2);
  col = 1:NbCol;
  col(end-1)=[];
  A = table2array(Tk(:,col)); %Remove the time stamp
  gof = A(:,end);
  [bestGOF, wBestGOF]=min(gof); %Find the minimum gof. And start from there
  k0 = A(wBestGOF(end),1:end-1) %Initial guess
  %k0 = Tk{end,1:end-1}
end

ub = ones (size(k0)).*200;
lb = zeros(size(k0));

%Optimise the rate constants
%----------------------------

%options = optimset('Display','iter','MaxFunEvals',100000, 'MaxIter',5000 , 'TolX' , 1e-12 , 'TolFun', 1e-12 );
%options = optimset('Display','iter', 'PlotFcns',@optimplotfval,'MaxFunEvals',100000, 'MaxIter',5000 ); %Fminsearch
options = optimset('Display','iter', 'PlotFcns', @optimplotfval , 'MaxFunEvals',100000, 'MaxIter',5000 ); %Simulated anealing


xflg = zeros(size(k0)); %Flag constant parameter
%optCase = 'O2';
optCase = 'Model';
%optCase = 'sigmoid';
%optCase = 'all';

switch optCase
  case 'Model'
    %fix the sigmoid parameter and optimise heuristic model
    %xflg(1) = 1; %Fix the G-factor
    xflg(9) = 1; %Fix the sigmoid gamma
    xflg(10) = 1; %Fix [LOOH]_50

    LOOHth = [];

  case 'sigmoid'
    %Fix heuristic model and optimise sigmoid shape
    xflg(1:8) = 1; %Fix the model parameters

  case 'all'
    %Vary all parameter, except the G-factor
    xflg(1) = 1; %Fix the G-factor

  case 'O2'
    xflg = ones(size(k0)); %Fix all
    xflg(6) = 0; %Except the reaction rate with O2

end

wvar = 1:numel(k0);
wvar(find(xflg)) = []; %Find the variable elements in the k vector
kin = k0(wvar) %the inital guess
ub = ub(wvar);
lb = lb(wvar);

%wp = [find(strcmp(Label , 'Sørensen et al., 2022'));find(strcmp(Label , 'Mascia et al., 2023'))];
%wp = find(strcmp(Label , 'Mascia et al., 2023'));
%wp = find(strcmp(Label , 'Tavakkoli et al., 2023'));


%wp = find((particle == 1) .* ~strcmp(Label , 'Duval et al., 2023')); %electron
%wp = find(~strcmp(Label , 'Duval et al., 2023'));

wp = find(~strcmp(Label , 'Duval et al., 2023') .* ~strcmp(Label , 'Konradsson et al., 2022') );
%wp = find(~strcmp(Label , 'Duval et al., 2023') .* ~strcmp(Label , 'Konradsson et al., 2022') .* ~strcmp(Label , 'Sesink et al. (2026)'));
%wp = find(strcmp(Label , 'Sesink et al. (2026)').* (numberOfPulses < 100)) ;
%wp = find(strcmp(Label , 'Sesink et al. (2026)')) ;
%wp = find((particle == 0) + strcmp(Label , 'Sesink et al. (2026)') ); %proton only + Sesink et al. (2026)

LabelList = unique({Label{wp}});

switch optCase
    case 'Model'
        LOOHth = [];
    case {'sigmoid', 'all' , 'O2'}

      if 1
        %Recompute [LOOHf] at each iteration
        %[~ , LOOHth] = GOFsimple(kin, k0 , xflg , O2(wp) , {DRp{wp}} ,{t{wp}} , no_repainting(wp) , particle(wp) , totalDose(wp) , numberOfPulses(wp) , SweepTime(wp) , DwellTime(wp), P_desq(wp)' ,{Label{wp}}, LabelList);
        %save(fullfile(folderData,'LOOHth.mat'),'LOOHth')
        LOOHth = [];
      else
        %Pre-compute [LOOHf]
        dataLoad = load(fullfile(folderData,'LOOHth.mat'));
        LOOHth = dataLoad.LOOHth;
      end

end

%For simulated anealing

%Optimise the parameters
%optimizer = 'fminsearch';
%optimizer = 'simulannealbnd';
optimizer = 'displayGOF';

switch optimizer
  case 'fminsearch'
    fileID = fopen(fileName,'w');
    [k,fval,exitflag] = fminsearch(@GOFsimple , kin , options , k0 , xflg , O2(wp) , {DRp{wp}} ,{t{wp}} , no_repainting(wp) , particle(wp) , totalDose(wp) , numberOfPulses(wp) , SweepTime(wp) , DwellTime(wp), P_desq(wp)' ,{Label{wp}}, LabelList,fileID,LOOHth);
    fclose(fileID);
    fprintf('gof = %3.8g \n' , fval)

  case 'simulannealbnd'
    %Include fixed parameters in function definition
    fileID = fopen(fileName,'w');
    Fobj = @(x)GOFsimple(x, k0 , xflg , O2(wp) , {DRp{wp}} ,{t{wp}} , no_repainting(wp) , particle(wp) , totalDose(wp) , numberOfPulses(wp) , SweepTime(wp) , DwellTime(wp), P_desq(wp)' ,{Label{wp}}, LabelList,fileID,LOOHth);
    [k,fval,exitflag] = simulannealbnd(Fobj , kin , lb , ub , options ); %Simulated enealing
    fclose(fileID);
    fprintf('gof = %3.8g \n' , fval)

  case 'displayGOF'
    %Display the graph with the latest value of the rate constant
    [gof , LOOHf] = GOFsimple(kin, k0 , xflg , O2(wp) , {DRp{wp}} ,{t{wp}} , no_repainting(wp) , particle(wp) , totalDose(wp) , numberOfPulses(wp) , SweepTime(wp) , DwellTime(wp), P_desq(wp)' ,{Label{wp}}, LabelList);
    fprintf('gof = %3.8g \n' , gof)
end

figure(10)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.021    0.122    0.45    0.77]);
saveas(10,fullfile(folderGraph,'complication_vs_LOOHf_sigmoid'),'jpg')
saveas(10,fullfile(folderGraph,'complication_vs_LOOHf_sigmoid'),'fig')
% saveas(2,fullfile(folderGraph,'complication_vs_LOOHf_sigmoid_zoom'),'jpg')


%----------------------------------------
% GOF to fit the simple function
%----------------------------------------
function [gof , LOOHf] = GOFsimple(pin, xcst , xflg , O2 , DRp ,t , no_repainting , particle , totalDose , numberOfPulses , SweepTime , DwellTime, Meas , Label, LabelList,fileID , LOOHth)

  if nargin >= 16
    logfile = true;
  else
    logfile = false;
  end

  k = buildVec(pin, xcst , xflg);
  k = abs(k);

  if logfile
      for i = 1:numel(k)-2
        fprintf(fileID,'%3.5g , ',k(i));
      end
  end


  %if prod(xflg(1:8)) & nargin >= 17 & ~isemty(LOOHth)
  if nargin >= 17 & ~isempty(LOOHth)
    %The model parameter are fixed
    %Use pre-computed [LOOH]f
    LOOHf = LOOHth;
  else
    %Compute the [LOOH]f
    LOOHf = getLOOHfMulti(k(1:8), O2 , DRp ,t , no_repainting , particle, totalDose , numberOfPulses , SweepTime , DwellTime);
  end

  gofType = 'sigmoid';
  %gofType = 'linear'

  switch gofType
    case 'sigmoid'
        %---- SIGMOID FIT -----
        %Compute the complication probability
        %Optimise the sigmoid parameters to fit the current data
        options = optimset('Display','off');
        ksig = fminsearch(@GOFSigmoid , [mean(LOOHf),mean(LOOHf)].*1e6 , options , LOOHf , Meas);
        ksig = abs(ksig); %Only accpt positive sigmoid parameters
        if logfile
            for i = 1:numel(ksig)
              fprintf(fileID,'%3.5g , ',ksig(i));
            end
        end

        param = abs([ksig(2).*1e-6,ksig(1).*1e-6,1,0]);
        P = sigmoid(param , LOOHf);
        RMS = sqrt(sum((P-Meas).^2)./numel(Meas));
        gof = RMS;

        xm = min(LOOHf):(max(LOOHf)-min(LOOHf))./100:max(LOOHf);
        xm = xm; %M
        ym = sigmoid(param , xm);


    case 'linear'
        % --- LINEAR FIT -----
        [b , xm , ym , r] = LinearReg(LOOHf' , Meas');
        gof = 1-abs(r(1,2));
    end

  if logfile
      fprintf(fileID,'%s , %3.8g \n',datetime('now'),gof);
  end

  %colour = colormap;
  symbols1 = {'o', 'square'	, 'diamond', 'pentagram', 'hexagram' 	, '^', 'v'	, '>'	, '<'	};
  symbols2 = {'+', 'x', '_', '|'	,  '.'};
  %colour={'r','y','g','c','b','m','k'};
  colour={'r','g','c','b','m','k'};


  figure(10)
  legendSTR = {};
  hold off
  S1 = 1;
  S2 = 1;
  for Sidx = 1:numel(LabelList)
    wType = find(strcmp(Label , LabelList{Sidx}));
    %colorN = getColIdx(DRp(wType) , min(DRp(wType)) , max(DRp(wType)))
    %scatter(LOOHf(wType)./LOOHf_50(wType) , Meas(wType) , [] , no_repainting(wType) , symbols{Sidx})
    colorIndex = no_repainting(wType) .* (no_repainting(wType) < 10) + 4 .* (no_repainting(wType) >= 10);
    O2Size = 2.*(1 + 100 .* (O2(wType) ./ (2 .* 50e-6))); %Size define oxygen level
    if particle(wType)
      %1= electron
      %open symbol
      scatter(LOOHf(wType).* 1e6 , Meas(wType) , O2Size , colorIndex , symbols2{S2})
      text(LOOHf(wType).* 1e6 , Meas(wType) ,[num2str(round(totalDose(wType))) repmat('Gy - ',size(totalDose(wType))) num2str(round(numberOfPulses(wType))) repmat('pls',size(numberOfPulses(wType)))])
      S2=S2+1;
    else
      %0= proton
      %close symbol
      scatter(LOOHf(wType).* 1e6 , Meas(wType) , O2Size , colorIndex , symbols1{S1},'filled')
      text(LOOHf(wType).* 1e6 , Meas(wType) ,[num2str(round(totalDose(wType))) repmat('Gy - ',size(totalDose(wType))) num2str(round(numberOfPulses(wType))) repmat('pls',size(numberOfPulses(wType)))])
      S1=S1+1;
    end
    hold on
    legendSTR{end+1} = Label{wType(1)};
    %text(LOOHf(wType(1))./LOOHf_50(wType(1)) , Meas(wType(1)),Label{wType(1)})
  end
  plot(xm .* 1e6 , ym,'-r')
  grid minor
  xlabel('[LOOH]_f = f(\int DR(t).dt , [O_2]_0) (MODEL)')
  ylabel('Complication (EXPERIMENT)')
  legend(legendSTR , 'location' , 'southeast')
  title('Colour: Nb paintings / Size : [O_2]')
  set(gca,'FontSize',20)
  drawnow


  % figure(2)
  % hold off
  % for Sidx = 1:numel(LabelList)
  %   wType = find(strcmp(Label , LabelList{Sidx}));
  %   %colorN = getColIdx(DRp(wType) , min(DRp(wType)) , max(DRp(wType)))
  %   scatter(LOOHf(wType)./LOOHf_50(wType) , Meas(wType) , symbols{Sidx})
  %   hold on
  %   %text(LOOHf(wType(1))./LOOHf_50(wType(1)) , Meas(wType(1)),Label{wType(1)})
  % end
  % plot(xm , ym,'-r')
  % grid minor
  % xlim([k(8).*0.5,k(8).*1.5])
  % xlabel('[LOOH]_f / [LOOH]_{LD50}  MODEL')
  % ylabel('Complication EXPERIMENT')
  % legend(legendSTR , 'location' , 'southeast')
  % drawnow

end

function gof = GOFSigmoid(k , LOOHf,Meas)
  param = abs([k(2).*1e-6,k(1).*1e-6,1,0]);
  P = sigmoid(param , LOOHf);
  RMS = sqrt(sum((P-Meas).^2)./numel(Meas));
  gof = RMS;
end


%---------------------------------
function p = buildVec(pin, xcst , xflg)
  p = xcst;
  varIdx = find(~xflg);
  p(varIdx) = pin;
end


%-------------------------------------
function LOOHf = getLOOHfMulti(k0, O2 , DRp ,t , no_repainting , particle, totalDose , numberOfPulses , SweepTime , DwellTime)
    %Compute the LOOHf
    ts1 = tic;
    for idx = 1:numel(particle)

      ts = tic;
      switch particle(idx)
        case 0 %PROTON
            %ts = tic;
            LOOHf(idx)  = getLOOHf2(k0, O2(idx) , DRp{idx} ,t{idx} , no_repainting(idx) , false);
            %fprintf('%d : [LOOH]f = %3.2g M T = %3.2g s \n',idx, LOOHf(idx),toc(ts))

        case 1 %PHOTON
          %LOOHf(idx)  = getLOOHf3(k0, O2(idx) , totalDose(idx) , numberOfPulses(idx) , SweepTime(idx) , DwellTime(idx), false);
          LOOHf(idx)  = getLOOHf4(k0, O2(idx) , totalDose(idx) , numberOfPulses(idx) , SweepTime(idx) , DwellTime(idx), false);
      end

      %fprintf('%d : %3.2g s : [LOOH]f = %3.2g M\n',idx,toc(ts),LOOHf(idx))
    end

end
