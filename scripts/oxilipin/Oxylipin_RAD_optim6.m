% Kinetic of creation of Oxylipin with radiation
% Optimise the endogenous rate constants to obtain the [Oxilipin] of Portier et al.
% Portier, L., Daira, P., Fourmaux, B. & Heinrich, S. Differential Remodeling of the Oxylipin Pool After FLASH Versus Conventional Dose-Rate Irradiation In Vitro and In Vivo. Int. J. Radiat. Oncol., Biol., Phys. 119, 1481–1492 (2024).
clear
close all

folder = 'D:\programs\openREGGUI\REGGUI_userdata\radiokinetics'

symbols = {'o', '*', '+',  '.', 'x', '_', '|', 'square'	, 'diamond'	, '^'	, 'v'	, '>'	, '<'	, 'pentagram', 'hexagram'};
lineType = { ':' , '-' , '--' , '-.'};
colors = {'k','b','g','r','c','m','y'};

%Collect all FLASH data into a single table
folder = 'C:\Users\rla\OneDrive - IBA Group\Projects\Flash\documents\biology\Temperature effect\experiments\Data_30_8_23'
fileName = 'Total Oxylipins - Temperature.xls'


DataSets{1} = {'RPE1FLASH' , 'RPE1CONV' , 'RPE1NI'};
DataSets{2} = {'SK28FLASH' , 'SK28CONV' , 'SK28NI'};
Celltype = {'RPE1-hTERT' , 'SK-MEL-28'};

dataXLS = readtable(fullfile(folder,fileName));
T = dataXLS(:,'Temperature').Variables;

for dIDX = 2:2

            [FLoxylipin , FLoxylipinSD] = getOxylip(dataXLS , DataSets{dIDX}{1});
            [CONVoxylipin , CONVoxylipinSD] = getOxylip(dataXLS ,DataSets{dIDX}{2});
            [NIoxylipin , NIoxylipinSD] = getOxylip(dataXLS , DataSets{dIDX}{3});

            %                        CONV                      UHDR
            data     =   [CONVoxylipin                   ; FLoxylipin];
            data_SD  =   [CONVoxylipinSD                 ; FLoxylipinSD];
            NbPulses =   [ones(size(CONVoxylipin)).* 250 ; ones(size(FLoxylipin)).* 5 ];
            Period = 1./ [ones(size(CONVoxylipin)).* 10  ; ones(size(FLoxylipin)).* 250 ] ; %Period = 1/ frequency
            PulseWidth = [ones(size(CONVoxylipin)).* 4   ; ones(size(FLoxylipin)).* 3.5 ] .* 1e-6;  %s
            O2         = ones(size(data)) .* 200;  %u-mol/l
            TotalDose =  ones(size(data)) .* 20; % Gy
            Tdata = [T ; T]; %°C

            figure(299+dIDX)
            legendSTR = {};
            errorbar(T, NIoxylipin , NIoxylipinSD, '^g')
            hold on
            legendSTR{end+1} = 'Non-irradiated';
            NIOxyL_r  = getOxiLBase(T, T , NIoxylipin)
            plot(T, NIOxyL_r , '-g' ) %non irradiated baseline
            legendSTR{end+1} = 'Non-irradiated';
            % legendSTR{end+1} = 'Non-irradiated';
            hold on


            %Plot the normalised experimental data
            figure(299+dIDX)
            hold on
            errorbar(T, CONVoxylipin , CONVoxylipinSD, 'ob')
            legendSTR{end+1} = 'CONV-RT';
            errorbar(T, FLoxylipin , FLoxylipinSD, '*r')
            legendSTR{end+1} = 'FLASH-RT';

            grid minor
            xlabel('Temperature (°C)')
            ylabel('Total [Oxylipin] (pmol per million cell)')
            title(Celltype{dIDX})


            %Get the basal level of oxylipin as a function of temperature
            [NIOxyL_r , slope , origin]  = getOxiLBase(Tdata, T , NIoxylipin);
            fprintf('Baseline : [OxyL]NI = %f + %f T \n', origin, slope );

            %Optimise the rate constants
            options = optimset('Display','iter','PlotFcns','optimplotfval');
            A = [];
            bs = [];
            Aeq = [];
            beq = [];
            nonlcon = [];


            %Load result for healthy cells as the starting point
            dataL = load(fullfile(folder,['oxilipin_optimum_k_' Celltype{1} '_C2.mat']));
            x0 = dataL.x;

                    %(kJ/mol)
                      %Ea        kbr2         kb3    Ea       kROOself    Ea        kb11       kbr     kb8    alpha
                      %   L*+L*           L* + O2              LOO*+LOO*           LOO*+LH
            isEa  =  [  1         0         0         1         0         1        0            0        0      0   ];
            Factor = [ 10    , 1e6       , 1e6     , 10     , 1e4     , 10        , 1      , 100    ,  1      , 1 ];
            lb     = [1e-1   , 1e-1      , 1       , 1e-1   , 1       , 1e-1      , 0.01   , 3e-1   ,  0.03 , 1e-2];  %Lower bound
            ub     = [10     , 1e3       , 3e2     , 10     , 3e3     , 10        , 100    , 4      ,  10   , 100] ; %Upper bound

            isEa = logical(isEa);
            lb(isEa) = x0(isEa); %Ea can only increase
            ub(isEa) = 20; %Ea is allowed to double
            ub(~isEa) = x0(~isEa); %k can only decrease
            ub(end) = 100;


            if (0)
              %Optimise rate constant
              [x,fval,exitflag]   = fmincon   (@GOFLOOHf,x0 ,A,bs,Aeq,beq,lb,ub,nonlcon , options , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata , data , Factor , NIOxyL_r);
              %[x,fval,exitflag]  = fminsearch(@GOFLOOHf,x0 ,options , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata , data , Factor , NIOxyL_r);

              x

              %save the results
              %save(fullfile(folder,['oxilipin_optimum_k_' Celltype{dIDX} '.mat']), '-v7.3')
              save(fullfile(folder,['oxilipin_optimum_k_' Celltype{dIDX} '_C2.mat']), '-v7.3')

            else
               %dataL = load(fullfile(folder,['oxilipin_optimum_k_' Celltype{dIDX} '.mat']));
               dataL = load(fullfile(folder,['oxilipin_optimum_k_' Celltype{dIDX} '_C2.mat']));

               x = dataL.x;

            end

            xi = abs(x .* Factor);
            fprintf('kb11 (25°C)     = %f \tEa = %f kJ/mol \n', xi(7) , xi(6))
            fprintf('kb3 (25°C)      = %1.2e \n', xi(3))
            fprintf('kbr2 (25°C)     = %1.2e \tEa = %f kJ/mol \n',xi(2) , xi(1) )
            fprintf('kROOself (25°C) = %1.2e \tEa = %f kJ/mol \n', xi(5) , xi(4))
            fprintf('kb8 (25°C)      = %f   \n', xi(9))
            fprintf('kbr (25°C)      = %f   \n', xi(8))
            fprintf('alpha           = %f   \n', xi(10))


            %Temperature dependence of rate Constants
            % High oxygen concentration
            %-------------------------------------------
            data2 = zeros(size(Tdata));

            [~ , OxyL] = GOFLOOHf(x , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata , data2 , Factor , NIOxyL_r , false);

            figure(299+dIDX)
            plot(Tdata(1:numel(T)), OxyL(1:numel(T)),'-b')
            legendSTR{end+1} = 'CONV  Model';
            plot(Tdata(numel(T)+1:numel(Tdata)), OxyL(numel(T)+1:numel(Tdata)),':r')
            legendSTR{end+1} = 'FLASH Model';
            legend(legendSTR)

end %for dIDX



%------------------------------------------------------
% The first data point must be the reference point : conventional DR at low temperature
function [gof , OxyL] = GOFLOOHf(x , TotalDose , Period , PulseWidth , NbPulses , O2 , T , data , Factor , NIOxyL_r  , flag)

      if nargin < 12
        flag = true;
      end

      x = abs(x)
      x = x .* Factor;
      alpha = x(10);

      %Get the predicted [LOOH]f
      for idx = 1:numel(TotalDose)
        kValue = getKt(T(idx) , x);
                     %getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue , verbose)
        LOOHf(idx)  = getLOOHf(TotalDose(idx) , Period(idx) , PulseWidth(idx) , NbPulses(idx) , O2(idx), kValue , false);
      end

      OxyL = alpha .* LOOHf + NIOxyL_r; %Add base line and multiply by alpha factor to obtain oxylipin concentration

      gof =  sum((data - OxyL').^2);

      if flag
        plotFit(T, data, OxyL);
      end

end

%------------------------------------------------------
function kValue = getKt(T , x)

  [kbr2 , kROOself , kb3 , kb8 , kbr] = getRateConstants();

  %diffusion in lipid membrane: activation energy
  kbr2     = k_arrhenius(T , x(2) , x(1).*1e3); % L* + L*
  kROOself = k_arrhenius(T , x(5) , x(4).*1e3); % ROO* + ROO*
  kb11     = k_arrhenius(T , x(7) , x(6).*1e3); % LOO* + LH

  %diffusion in liquid: diffusion constant
  kb8  = kDT(x(9),T); %LOO* + XSH -> LOOH + XS*
  kbr  = kDT(x(8),T); %L*  + XSH -> LH v+ XS*
  kb3  = kDT(x(3),T); %R* + O2

  kValue = [ kbr2 ,    kb3 ,  kROOself  , kb8 , kbr , kb11];

end

%------------------------------------------------------
% Variation of rate constant due to diffusion
%------------------------------------------------------
function k = kDT(k,T)
  k = k .* (T + 273) ./ (25 + 273);
end

%-------------------------------------------------
% Compute the normalised [LOOH]f expected from basal metabolism
% whithout irradiation
% as a function of temperature
%
% NIOxyL_r : Oxilipin concentration at that temperature. Normalised to 1 at minimum temeprature
%-------------------------------------------------
function [NIOxyL_r , slope , origin]  = getOxiLBase(Tinter, T , NIoxylipin)

  %Basal level of oxylipin
  NIOxyL_r = NIoxylipin(1) + (Tinter-T(1)) .* (NIoxylipin(4)-NIoxylipin(1)) ./ (T(4)-T(1));
  NIOxyL_r = NIOxyL_r';
  origin = NIoxylipin(1) - T(1) .* (NIoxylipin(4)-NIoxylipin(1)) ./ (T(4)-T(1));
  slope = (NIoxylipin(4)-NIoxylipin(1)) ./ (T(4)-T(1));

end


%---------------------------------------------------------------------------------
function [oxylipin , oxylipinSD] = getOxylip(data , ColmunName)


    oxylipin = data(:,ColmunName).Variables; %Oxylipin concentration

    colNb = find(strcmp(data.Properties.VariableNames , ColmunName)) + 1;
    oxylipinSD = data(:,data.Properties.VariableNames{colNb}).Variables; %standard deviation on concentration


end

%---------------------------------------------------------------
% Compute the ratio (A-B)/(C-D)
% Compuyte the error bar for the function F= (A-B)/(C-D), knowing error bar on A,B,C,D
%
% error = sqrt( sum_u (dF/du .* Su).^2  )
%
% INPUT
%---------------------------------------------------------------
function [F , SD_F] = getRatio(A , B , C , D , sA , sB , sC , sD)

  F = (A-B) ./ (C-D); %The ratio
  SD_F = (sA .* 1 ./ (C-D)).^2 + (-sB .* 1 ./ (C-D)).^2 + (sC .* (A-B)./(C-D).^2).^2 + (-sD .*(A-B)./(C-D).^2).^2;
  SD_F = sqrt(SD_F); %Error bars computed using formula for propagation of uncertainty
  %https://en.wikipedia.org/wiki/Propagation_of_uncertainty

end


function plotFit(Tdata, data, model)
  figure(400)
  hold off
  plot(Tdata(1:4), data(1:4),'ob')
  hold on
  plot(Tdata(1:4), model(1:4),'-b')
  plot(Tdata(5:8), data(5:8),'or')
  plot(Tdata(5:8), model(5:8),'-r')
  grid minor
  xlabel('Temperature (°C)')
  ylabel('Total [Oxylipin] (pmol per million cell)')
  drawnow
end
