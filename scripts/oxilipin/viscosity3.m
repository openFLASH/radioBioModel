%Effect of temperature on viscosity


%%REFERENCE
% [1] https://en.wikipedia.org/wiki/Temperature_dependence_of_viscosity
% [2] Bag, N., Hui, D., Yap, X. & Wohland, T. Temperature dependence of diffusion in model and live cell membranes characterized by imaging fl uorescence correlation spectroscopy. Biochim. Biophys. Acta 1838, 802–813 (2014).
% [3] DATA BOOK ON THE VISCOSITY OF LIQUIDS D.s. Viswanath, Author ; G. Natarajan	Hemisphere Publishing Corporation, 1989 ISBN 978-0-89116-778-5

close all
clear

folderGraph ="C:\Users\rla\IBA Group\Proton Flash Therapy - Documents\Local effect model to describe FLASH\temperature_viscosity\figures\data";
fileName ="C:\Users\rla\IBA Group\Proton Flash Therapy - Documents\Local effect model to describe FLASH\temperature_viscosity\data\k.txt";

symbols = {'o', '*', '+',  '.', 'x', '_', '|', 'square'	, 'diamond'	, '^'	, 'v'	, '>'	, '<'	, 'pentagram', 'hexagram'};
lineType = { ':' , '-' , '--' , '-.'};
colors = {'k','b','g','r','c','m','y'};


varName = {'kdr'    , 'Ea'  , 'kb2DR' , 'kb3_298' , 'kscale'};
%k     = [  0.010092 ,  52   , 77.543  , 0.010275 , 99.313];
k     = [  0.1        ,  52   ,  1  , 0.1 , 1];

%Combined effect of activation energy + viscosity D298K
%----------------------------
T = 25:40;
T = T + 273;

kdr = k(1);

%Ea = 50 .* ones(1,5);
%k298K = [1e-2 1e-1 1 10 100];

%Ea = [10 , 36 , 37, 50];
%k298K = ones(size(Ea));

Ea  = 10:5:50; %kJ/mol
k298K = getVisMem(Ea).*k(5); %Use cell data to predict viscosity from activation energy

kb2DR = k(3);
kb3_298 = k(4);


for idx = 1:numel(Ea)

      [LOOr , kb3 , kd2r] = getLOORatio(kdr , T , kb3_298 , Ea(idx) , k298K(idx) , kb2DR);

      figure(17)
      plot(T,LOOr)
      text(T(end) , LOOr(end) ,['Ea=' num2str(Ea(idx)) , ' k_{289k} = ' , num2str(k298K(idx)) ])
      text(T(1) , LOOr(1) ,['Ea=' num2str(Ea(idx)) , ' k_{289k} = ' , num2str(k298K(idx)) ])
      hold on

      if idx == 1
          figure(18)
          plot(T,kb3)
          text(T(end) , kb3(end) ,['Ea=' num2str(Ea(idx)) , ' k_{289k} = ' , num2str(k298K(idx)) ])
          hold on

          figure(19)
          plot(T,kb3 ./ min(kb3) , '--r')
          text(T(end) , kb3(end)./kb3(1) ,'Water')
          % plot(T,kb3 )
          % text(T(end) , kb3(end) ,'Water')
          hold on
      end

      figure(19)
      plot(T,kd2r ./ min(kd2r))
      text(T(end) , kd2r(end)./ min(kd2r) ,['Ea=' num2str(Ea(idx)) , ' k_{289k} = ' , num2str(k298K(idx)) ])
      %plot(T,kd2r )
      %text(T(end) , kd2r(end) ,['Ea=' num2str(Ea(idx)) , ' k_{289k} = ' , num2str(k298K(idx)) ])

      hold on

end

figure(17)
grid minor
xlabel('T (K)')
ylabel('[LOO^*]_r')
title('E_a and D_{298K}')

figure(18)
grid minor
xlabel('T (K)')
ylabel('k_{b3}')
title('k_{b3}')

figure(19)
grid minor
xlabel('T (K)')
ylabel('k_{d2r}')
title('k_{d2r}')


%-----------------------------------
% Goodness of fit function
% k = [kdr , Ea , k298K , kb2DR , kb3_298]
%-----------------------------------
function gof = gofLOOR(kin, k0, xflg , fileID)

  k = buildVec(kin, k0 , xflg);
  k = abs(k);
  kdr = k(1);
  Ea = k(2);
  k298K = k(3);
  kb2DR = k(4);
  kb3_298 = k(5);

  k298K = k(3);

  T = [25,38] + 273;

  LOOr = getLOORatio(kdr , T , kb3_298 , Ea , k298K , kb2DR);


  gof = LOOr(2) - LOOr(1); %We want this to be as negative as possible

  for i = 1:numel(k)
    fprintf(fileID,'%3.5g , ',k(i));
  end
  fprintf(fileID,'%3.8g \n',gof);

end

function gof = gofLOOR2(kin, k0, xflg , fileID)

  k = buildVec(kin, k0 , xflg);
  k = abs(k);

  kdr = k(1);
  Ea = k(2);
  kb2DR = k(3);
  kb3_298 = k(4);
  k298K = getVisMem(Ea) .* k(5);

  T = [25,38] + 273;

  LOOr = getLOORatio(kdr , T , kb3_298 , Ea , k298K , kb2DR);


  gof = LOOr(2) - LOOr(1); %We want this to be as negative as possible

  for i = 1:numel(k)
    fprintf(fileID,'%3.5g , ',k(i));
  end
  fprintf(fileID,'%3.8g \n',gof);

end

%-------------------------------------------
% Compute [LOOR]
%-------------------------------------------
function LOOR = getLOOR(kdr,kb3,kd2r,kb2DR)
  F3 = ThirdFactor(kdr,kb3,kd2r,kb2DR);
  F2 = SecondFactor(kb3,kd2r);
  LOOR = F2 .* F3;
end

% ------------------------------------
% Compute the third factor in the product
% to ocmpute the [LOOr] at steady state
% ------------------------------------
function F3 = ThirdFactor(kdr,kb3,kd2r,kb2DR)
  F3 = sqrt( (kdr + kb3).^2 + 4.* kd2r.* kb2DR ) - (kdr + kb3);
end

% ------------------------------------
% Compute the second factor in the product
% to ocmpute the [LOOr] at steady state
% ------------------------------------
function F2 = SecondFactor(kb3,kd2r)
  F2 = kb3 ./ kd2r;
end

%-----------------------------------
% kb3 as a function of temperature
%-----------------------------------
function kb3 = g_kb3_T(T, kb3_298)

  nu298 = viscTempE(298);
  nu    = viscTempE(T);

  kb3 = kb3_298 .* T .* nu298 ./ (nu .* 298);
  %kb3 = kb3_298 .* T  ./ 298;


end

%-----------------------------------
% kd2r as a function of temperature
%-----------------------------------
function kd2r = g_kd2r_T(T, Ea , k298K)
  kd2r = Dlip(T, Ea , k298K);
end


%------------------------------
% Empirical description of viscosity of water vs temperature
% [1]
%------------------------------
function nu = viscTempE(T)
  %Water [1]
  cofW(1,1) = 1.856e-11; %mPa.s
  cofW(1,2) = 4209; %K
  cofW(1,3)= 0.04527; %K-1
  cofW(1,4)= -3.376e-5; %K-2
  nu = cofW(1) .* exp(cofW(2) ./ T + cofW(3) .* T + cofW(4) .* T.^2); %Viscosity (mPa.s)
end

%---------------------------------
%Diffusion coefficient in lipd
%---------------------------------
function D = Dlip(T, Ea , D298K)
  R = 8.31446261815324; % J mol−1 K−1
  D = D298K .* exp(-Ea.*1e3.*(1./T-1./298)./R);
end

%---------------------------------
function p = buildVec(pin, xcst , xflg)
  p = xcst;
  varIdx = find(~xflg);
  p(varIdx) = pin;
end

%-------------------------------
% Compute the [LOO*] ratio
%-------------------------------
function [LOOr , kb3 , kd2r] = getLOORatio(kdr , T , kb3_298 , Ea , k298K , kb2DR)

      kb3 =  g_kb3_T(T, kb3_298);
      kd2r = g_kd2r_T(T, Ea , k298K);

      LOOr = getLOOR(kdr,kb3,kd2r,kb2DR);
      LOOr = LOOr ./ max(LOOr);

end

%---------------------------------------------
% Get approxiamte properties of membranes
%---------------------------------------------
function k298K = getVisMem(Ea)
  DEacoef = [-0.0583 , 3.2]';
  k298K = DEacoef' * [Ea ; ones(size(Ea))];
end
