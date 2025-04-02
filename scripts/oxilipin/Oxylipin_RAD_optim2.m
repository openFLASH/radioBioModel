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

data = readtable(fullfile(folder,fileName));
T = data(:,'Temperature').Variables;
[FLoxylipin , FLoxylipinSD] = getOxylip(data , 'RPE1FLASH');
[CONVoxylipin , CONVoxylipinSD] = getOxylip(data , 'RPE1CONV');
[NIoxylipin , NIoxylipinSD] = getOxylip(data , 'RPE1NI');

%                        CONV                      UHDR
data     =   [CONVoxylipin                   ; FLoxylipin];
data_SD  =   [CONVoxylipinSD                 ; FLoxylipinSD]
NbPulses =   [ones(size(CONVoxylipin)).* 250 ; ones(size(FLoxylipin)).* 5 ];
Period = 1./ [ones(size(CONVoxylipin)).* 10  ; ones(size(FLoxylipin)).* 250 ] ; %Period = 1/ frequency
PulseWidth = [ones(size(CONVoxylipin)).* 4   ; ones(size(FLoxylipin)).* 3.5 ] .* 1e-6;  %s
O2         = ones(size(data)) .* 200;  %u-mol/l
TotalDose =  ones(size(data)) .* 20; % Gy
Tdata = [T ; T]; %°C

%Get the value to renormalise
% [RefOx , RefOx_idx] = max(fullArray,[],'all') %Get the value ofthe largest concentration
% RefOxSD = fullArraySD(RefOx_idx) %Get the SD on the largest concentration
[~ , Tmin_idx] = min(T);
RefOx = CONVoxylipin(Tmin_idx);
RefOxSD = CONVoxylipinSD(Tmin_idx);

%Plot the normalised experimental data
figure(300)
hold on
errorbar(T, FLoxylipin , FLoxylipinSD, '*r')
errorbar(T, CONVoxylipin , CONVoxylipinSD, 'ob')
errorbar(T, NIoxylipin , NIoxylipinSD, '^g')

grid minor
xlabel('Temperature (°C)')
ylabel('[Oxylipin] nmol/g')

figure(300)
legendSTR = {};
NIOxyL_r  = getOxiLBase(T, T , NIoxylipin)
plot(T, NIOxyL_r , ['-' colors{3}]) %FLASH
legendSTR{end+1} = 'Non irradiated';
hold on


%Get the basal level of oxylipin as a function of temperature
[NIOxyL_r , slope , origin]  = getOxiLBase(Tdata, T , NIoxylipin);

%Normalize with the radio induced [] at CONV and 20°C
[data , SD_data] = getRatio(data , NIOxyL_r' , data(1) , NIOxyL_r(1) , data_SD , zeros(size(NIOxyL_r')) , data_SD(1) , 0)

%Compute the partition functions
% Hp  = power(10,-7); %[H+] Buffered solution -> constant concentration
% [PF.OHr , PF.Orm, PF.dOHr_dCt , PF.dOrm_dCt] = acidPartition(1 , Hp , 11.9);
% [PF.H2O2 , PF.HO2m, PF.dH2O2_dCt , PF.dHO2m_dCt] = acidPartition(1 , Hp , 11.7);
% [PF.HO2r , PF.O2rm, PF.dHO2r_dCt , PF.dO2rm_dCt] = acidPartition(1 , Hp , 4.9);



%Optimise the rate constants
options = optimset('Display','iter','PlotFcns','optimplotfval');
A = [];
bs = [];
Aeq = [];
beq = [];
nonlcon = [];
lb = [1     , 1e5  ,  1e2 , 1    , 1e4];  %Lower bound
ub = [1000  , 3e9  ,  3e8 , 1000 , 3e7] ; %Upper bound
        %Ea   kbr2    kb3    Ea   kROOself
%             R*+R*   R* + O2

        %(kJ/mol)
          %Ea   kbr2    kb3    Ea   kROOself    Ea    kb11  kbr     kb8
          %   L*+L*   L* + O2   LOO*+LOO*       LOO*+LH
%x0     = [ 19 , 21  , 46       , 3   , 2       , 1   , 20  , 1.25 ,  5]
x0     = [24.5122   23.4970   12.8920    7.7660    3.1590    1.2650   15.1925    0.7564    9.6434]
Factor = [ 10 , 1e6 , 1e6      , 100 , 1e4     , 100 , 20  , 100 ,  1e-2];

if (0)
  %Optimise rate constant
  %[x,fval,exitflag] = fmincon(@GOFLOOHf,x0 ,A,bs,Aeq,beq,lb,ub,nonlcon , options , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata , data , PF);
  [x,fval,exitflag] = fminsearch(@GOFLOOHf,x0 ,options , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata , data , Factor);

  x

  %save the results
  save(fullfile(folder,['oxilipin_optimum_k.mat']), '-v7.3')
else
   dataL = load(fullfile(folder,['oxilipin_optimum_k.mat']));
   x = dataL.x; %TODO !!!!
  %x= [25.2005   25.9819   13.9496    7.2889    3.4262    1.2185   17.8155    0.8029    6.4052];
end



%Temperature dependence of rate Constants
% High oxygen concentration
%-------------------------------------------
T = 20:37; %°C
T = T';
%                        CONV                      UHDR
NbPulses =   [ones(size(T)).* 250 ; ones(size(T)).* 5 ];
Period = 1./ [ones(size(T)).* 10  ; ones(size(T)).* 250 ] ; %Period = 1/ frequency
PulseWidth = [ones(size(T)).* 4   ; ones(size(T)).* 3.5 ] .* 1e-6;  %s
Tdata2 = [T ; T]; %°C
O2         = ones(size(Tdata2)) .* 200;  %u-mol/l
TotalDose =  ones(size(Tdata2)) .* 20; % Gy
data2 = zeros(size(Tdata2));

[~ , LOOHf] = GOFLOOHf(x , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata2 , data2 , Factor , false);
LOOHfr =  LOOHf ./ LOOHf(1); %Compute the relative oxylipin concentrations

figure(400)
legendSTR={};
hold off
errorbar(Tdata(1:4), data(1:4) , SD_data(1:4), 'ob')
legendSTR{end+1} = 'CONV  Measurement';
hold on
plot(Tdata2(1:numel(T)), LOOHfr(1:numel(T)),'-b')
legendSTR{end+1} = 'CONV  Model';
errorbar(Tdata(5:8), data(5:8) , SD_data(5:8), 'or')
legendSTR{end+1} = 'FLASH Measurement';
plot(Tdata2(numel(T)+1:numel(Tdata2)), LOOHfr(numel(T)+1:numel(Tdata2)),'-r')
legendSTR{end+1} = 'FLASH Model';

grid minor
xlabel('Temperature (°C)')
ylabel('[LOOHf]_f / [LOOHf]_f(conv,20°C)')

pause


% %Reducing diffusion rate in lipid membrane
% %-------------------------------------------
% %Ea   kbr2    kb3    Ea   kROOself    Ea    kb11  kbr     kb8
% %   L*+L*   L* + O2   LOO*+LOO*       LOO*+LH     L*+AH   LOO*+ AH
Mask_K = zeros(1,9);
KMemb = [2, 5 , 7];
Mask_K(KMemb) = 1; %Identify the rate constants that would be different in tumor cells
S_k = 1/2.6
x_slow = x .* Mask_K .* S_k + x .* ~Mask_K;
%
% [~ , LOOHf_s] = GOFLOOHf(x_slow , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata2 , data2 , Factor, false);
%
% LOOHf_s =  LOOHf_s ./ LOOHf(1); %Compute the relative oxylipin concentrations
%
% figure(400)
% hold on
% plot(Tdata2(1:numel(T)), LOOHf_s(1:numel(T)),'--b')
% plot(Tdata2(numel(T)+1:numel(Tdata2)), LOOHf_s(numel(T)+1:numel(Tdata2)),'--r')
% legendSTR{end+1} = 'CONV  Model with k/2.6';
% legendSTR{end+1} = 'FLASH Model with k/2.6';
% pause
%
% %Increasing the activation energy in memebrane
% %-------------------------------------------
% %Ea   kbr2    kb3    Ea   kROOself    Ea    kb11  kbr     kb8
% %   L*+L*   L* + O2   LOO*+LOO*       LOO*+LH     L*+AH   LOO*+ AH
Mask_E = zeros(1,9);
KMemb = [1, 4 , 6];
Mask_E(KMemb) = 1; %Identify the rate constants that would be different in tumor cells
S_E = 0.5
x_slow = x .* Mask_E .* S_E + x .* ~Mask_E;
%
% [~ , LOOHf_s2] = GOFLOOHf(x_slow , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata2 , data2 , Factor, false);
%
% LOOHf_s2 =  LOOHf_s2 ./ LOOHf(1); %Compute the relative oxylipin concentrations
%
% figure(400)
% hold on
% plot(Tdata2(1:numel(T)), LOOHf_s2(1:numel(T)),'-.b')
% plot(Tdata2(numel(T)+1:numel(Tdata2)), LOOHf_s2(numel(T)+1:numel(Tdata2)),'-.r')
% legendSTR{end+1} = 'CONV  Model with Ea/2';
% legendSTR{end+1} = 'FLASH Model with Ea/2';
% pause

%Increasing the activation energy & reducing diffusion in memebrane
%-------------------------------------------
%Ea   kbr2    kb3    Ea   kROOself    Ea    kb11  kbr     kb8
%   L*+L*   L* + O2   LOO*+LOO*       LOO*+LH     L*+AH   LOO*+ AH
x_slow = x .* Mask_E .* S_E + x .* Mask_K .* S_k + x .* ~(Mask_E+Mask_K);

[~ , LOOHf_s1] = GOFLOOHf(x_slow , TotalDose , Period , PulseWidth , NbPulses , O2 , Tdata2 , data2 , Factor, false);

LOOHf_s1 =  LOOHf_s1 ./ LOOHf(1); %Compute the relative oxylipin concentrations

figure(400)
hold on
plot(Tdata2(1:numel(T)), LOOHf_s1(1:numel(T)),':b','Linewidth', 2)
plot(Tdata2(numel(T)+1:numel(Tdata2)), LOOHf_s1(numel(T)+1:numel(Tdata2)),':r','Linewidth', 2)
legendSTR{end+1} = 'CONV  Model with k/2.6 & Ea/2';
legendSTR{end+1} = 'FLASH Model with k/2.6 & Ea/2';

legend(legendSTR , 'Location' , 'eastoutside')
pause

%------------------------------------------------------
% The first data point must be the reference point : conventional DR at low temperature
function [gof , LOOHf] = GOFLOOHf(x , TotalDose , Period , PulseWidth , NbPulses , O2 , T , data , Factor , flag)

      if nargin < 11
        flag = true;
      end

      x = abs(x)
      x = x .* Factor;

      for idx = 1:numel(TotalDose)
        kValue = getKt(T(idx) , x);
                     %getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2 , kValue , verbose)
        LOOHf(idx)  = getLOOHf(TotalDose(idx) , Period(idx) , PulseWidth(idx) , NbPulses(idx) , O2(idx), kValue , false);
      end

      LOOHfr = LOOHf ./ LOOHf(1); %Normalise by the first data point
      gof =  sum((data - LOOHfr').^2);

      if flag
        plotFit(T, data, LOOHfr);
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
  kb8  = kDT(x(8),T); %LOO* + XSH -> LOOH + XS*
  kbr  = kDT(x(9),T); %L*  + XSH -> LH v+ XS*
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
  ylabel('[LOOHf]_f / [LOOHf]_f(conv,20°C)')
  drawnow
end
