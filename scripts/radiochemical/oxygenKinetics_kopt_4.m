% Optimise the value of the rate constants that maximise the correlation between [LOOH]f and the crypt survival
% AND maintain the decrease of oxygen concentration at 0.45 umol/l/Gy
% Use the radio-kinetic model for the 2 phases model to compute the concentration of different radical species as a function of time
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

clear
close all

fileName = 'D:\programs\openREGGUI\REGGUI_userdata\radiokinetics\x.csv';

options = optimset('Display','iter');
    %[kbr2     , kb3  ,  kLOOself   , kb8     , kbr]
x0 = [1    1   1  1  1]; %Initial guess Multiplicative factor of the current rate constant

[kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants();
Kref = [kbr2       , kb3       ,  kLOOself       , kb8 , kbr];


A = [];
bs = [];
Aeq = [];
beq = [];
nonlcon = [];


%NB: contrarily to what is suggested on the MAtlab documentation, the boundary conditions are defined as STRICTLY larger or STRICTLY smaller
    %[kbr2  , kb3  ,  kLOOself   , kb8     , kbr]
lb = [1e5   , 1e6  ,  1e4        , 0.03    , 39  ]  ;  %Lower bound
ub = [1e9   , 3e8  ,  3e7        , 10      , 400 ,] ; %Upper bound

%kbr2: [1] R* recombine via a second order reaction with a rate constant 10^5 < kd2R < 10^9 (mol/l s)^-1 consistent with the rate of lateral diffusion of lipids in biological membranes.
%kb3: [1] The rate constant approaches diffusion-limit [11,6], ranging from k 10^8 to 10^10 (mol/l s)^-1 in water, depending on the electronic properties of R. However, in the cellular environment, the diffusion rate is lower due to higher viscosity of cytosol: 10^6 to 10^7 (mol/l s)^-1 [8,23]
%     [2a,3a,5a,6a]
% kLOOself : [2a,3a,5a,6a]
%kb8 :  [2a,3a,5a,6a]

%[1] Labarbe, R., Hotoiu, L., Barbier, J. & Favaudon, V. A physicochemical model of reaction kinetics supports peroxyl radical recombination as the main determinant of the FLASH effect. Radiother. Oncol. 153, 303–310 (2020).
%[11] Gray B, Carmichael AJ. Kinetics of superoxide scavenging by dismutase enzymes and manganese mimics determined by electron spin resonance. Biochem J 1992;281:795–802.
%[6] Buxton, G V, Greenstock, C L, Helman, W P, Ross, A B, Tsang, W. Critical review of rate constants for reactions of hydrated electrons, hydrogen atoms and hydroxyl radicals (?OH/O-) in aqueous solution. J Phys Chem Ref Data, 17:513– 886, 1988.
%[8] Epp ER, Weiss H, Ling CC. Irradiation of cells by single and double pulses of high intensity irradiation: oxygen sensitization and diffusion kinetics. Curr Top Radiat Res 1976;11:201–30.
%[23] Michael BD, Davies S, Held KD. Ultrafast chemical repair of DNA single and double strand break precursors in irradiated V79 cells. Basic Life Sci 1986;38:89–100.
% [2a]	Stark, G. The effect of ionizing radiation on lipid membranes. Biochim. Biophys. Acta 1071, 103–122 (1991).
% [3a].	Babbs, C., and Steiner, M. G. (1990) Simulation of free radical reactions in biology and medicine: a new two-compartment kinetic model of intracellular lipid peroxidation. Free Radical Biol. Med. 8, 471-485.
% [5a].	Salvador, A. et al. Kinetic Modelling of in Vitro Lipid Peroxidation Experiments - ’ Low Level ’ Validation of a Model of in Vivo Lipid Peroxidation KINETIC MODELLING OF IN VZTRO LIPID. Free Rad. Rex 23, 151–172 (1995)
% [6a].	Antunes, F., Salvador, A., Marinho, H. S., Alves, R. & Ruy E Pinto. Lipid peroxidation in mitochondrial inner membranes . I . An integrative kinetic model. Free Radic. Biol. Med. 21, 917–943 (1996).

lb = lb ./ Kref
ub = ub ./ Kref


fileID = fopen(fileName,'w');
[x,fval,exitflag] = fmincon(@GOFcorrelation,x0 ,A,bs,Aeq,beq,lb,ub,nonlcon , options , fileID);
fclose(fileID);

for id = 1:numel(x)
  fprintf('k = %f \n',x(id))
end

%=================================================
% Goodness of fit function
% correlation between [LOOH]f and the crypt survival
%=================================================
function Rsq = GOFcorrelation(x , fileID)

      x = abs(x) %Ratye constant are positive

      [kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants();
              %[   kbr2    , kb3       ,   kLOOself       , kb8        , kbr]
      kValue = [x(1).*kbr2 , x(2).*kb3 ,  x(3).*kLOOself  , x(4).*kb8 , x(5).*kbr];

      % crypt survival for 11.2 Gy
      %-----------------------------
      O2         = [50      ,  50    ,  50   , 50    ,  50   ,  50   , 50   ,  50    , 50    , 50]; %u-mol/l
      TotalDose =  [11.2    , 11.2   ,  11.2 , 11.2  ,  11.2 ,  11.2 , 11.2 , 11.2   , 11.2  , 11.2] ;% Gy
      NbPulses   = [1       ,   2    ,   2   ,  2    ,   2   ,   2   , 5     , 30    , 100   , 300];
      DRp        = [3.39e6 , 1.6e6  , 1.6e6 , 1.6e6 , 1.6e6 , 1.6e6 , 6.6e5 , 1.1e5 , 3.2e4 , 1.1e4 ]; % Gy/s Peak dose rate
      Dp         = [11.2    , 5.6    , 5.6   , 5.6   , 5.6   , 5.6   , 2.24  , 0.37  , 0.112 , 0.037]; %Gy Dose per pulse
      PulseWidth = Dp ./ DRp ; %s
      Period     = [3.4e-6  , 3.3e-3 , 10e-3 , 40e-3 ,   3   ,  30   , 10e-3 , 10e-3 , 10e-3 , 10e-3 ]; %Period (s)
      Label      = {};
      Meas       = [27.6    ,  14.9  , 13.2  , 6.9   , 8.1   , 6.4   , 13.8  , 11.5  , 6.7   , 8.6 ]; %Crypt survival Table 1
      SD         = [ 4      ,  2.2   , 1.6   , 2.8   , 2.5   , 1.5   , 0.8   , 2.2   , 1.8   , 1.3 ] ; %Error bars
      doseRateEXP= [3.3e6   , 3.4e3  , 1.1e3 , 280   , 3.7   , 0.37  , 280   , 39    , 11    , 3.7 ]; %Gy/s Average dose rate reported in paper
                    %NB: The average dose rate reported in the article does not seem correct.
      % crypt survival for 11.2 Gy whole abdominal irradiation of mice aged 9 to 10 weeks
      %Ruan, J. et al. Irradiation at Ultra-High ( FLASH ) Dose Rates Reduces Acute Normal Tissue Toxicity in the Mouse Gastrointestinal System. Int. J. Radiat. Oncol. Biol. Phys. 111, 1250–1261 (2021).


      %Compute the [LOOH]f for all the experimental conditions
      [LOOHf , AvDoseRate , PkDoseRate , O2f ] = getLOOHf(TotalDose , Period , PulseWidth , NbPulses , O2, kValue , false);


      % crypt survival for 12.5 Gy
      %-----------------------------
      O22         = [50      ,  50    ,  50   , 50    ,  50   ,  50   , 50    ,  50    , 50    , 50    , 50   , 50]; %u-mol/l
      TotalDose2 =  [12.5    , 12.5   ,  12.5 , 12.5  ,  12.5 ,  12.5 , 12.5  , 12.5   , 12.5  , 12.5  , 12.5 , 12.5] ;% Gy
      NbPulses2   = [1       ,   2    ,   2   ,  2    ,   2   ,   2   , 5     , 10     , 30    , 100   , 300  , 1250];
      DRp2        = [3.7e6   , 1.8e6  , 1.8e6 , 1.8e6 , 1.8e6 , 1.8e6 ,7.4e5 , 3.2e5   , 1.2e5 , 3.7e4 , 1.2e4, 2.9e3 ]; % Gy/s Peak dose rate
      Dp2         = [12.5    , 6.25    , 6.25 , 6.25  , 6.25  , 6.25  , 2.5  , 1.25    , 0.417 , 0.125 , 0.042, 0.01]; %Gy Dose per pulse
      PulseWidth2 = Dp2 ./ DRp2 ; %s
      Period2     = [3.4e-6  , 3.3e-3 , 10e-3 , 40e-3 ,   3   ,  30   , 10e-3 , 10e-3 , 10e-3 , 10e-3  , 10e-3, 40e-3 ]; %Period (s)
      Meas2       = [11      , 9.1    , 16.1  , 15.3  , 3.8   , 4.8   , 7.8   , 3.3   , 5.0   , 6.1    , 3.3  , 1.3 ]; %Crypt survival Table 1
      SD2         = [2.8     , 2.9    , 5.8   , 3.2   , 0.2   , 0.9   , 1.4   , 1.0   , 1.8   , 2.3    ,1.1   , 0.3 ] ; %Error bars
      % crypt survival for 12.5 Gy whole abdominal irradiation of mice aged 30 to 31 weeks
      %Ruan, J. et al. Irradiation at Ultra-High ( FLASH ) Dose Rates Reduces Acute Normal Tissue Toxicity in the Mouse Gastrointestinal System. Int. J. Radiat. Oncol. Biol. Phys. 111, 1250–1261 (2021).


      %Compute the [LOOH]f for all the experimental conditions
      [LOOHf2 , AvDoseRate2 , PkDoseRate2 , O2f2 ] = getLOOHf(TotalDose2 , Period2 , PulseWidth2 , NbPulses2 , O22, kValue , false);



      %Draw the regression line
      x = [LOOHf' ; LOOHf2'];
      x = [x , ones(size(x,1),1)];
      y = [Meas' ; Meas2'];
      b = x\y; %Make a linear regression of [LOOH]f vs log(DR)
      Predict = b'*x';

      %The GOF has one term that maximise the correclation between [LOOH]f and the biological output
      Rsq1 = -abs(1 - sum((y - Predict').^2) / sum((y - mean(y)).^2));
            % Coefficient of determination,
            %Put a negative number. fmincon search the minimum
            %but we do not want a correlation equal to zero. We want a correlation to be +1 or -1

      %The GOF has a second term that constraint the decrease of [O2] with a g-factor < 0.45 umol/l/Gy
      gO2 = ([O2 , O22] - [O2f , O2f2]) ./ [TotalDose , TotalDose2];   %u-mol/l / Gy G-factor for oxygen consumption
      Rsq2 = max( (gO2 - 0.45) , 0);  %we want the number to be <0.45 u-mol/l. put a cost when it is larger
                                     %Weiss H, Epp ER, Heslin JM, Ling CC, Santomasso A. Oxygen depletion in cells irradiated at ultra-high dose-rates and at conventional dose-rates. Int J Radiat Biol 1974;26:17–29.
      Rsq2 = mean(Rsq2);    %The GOF is the mean error, so it does not depend on the number of data points

      Rsq = Rsq1 + Rsq2;

      %Save current X to file
      for i = 1:numel(kValue)
        fprintf(fileID,'%3.5f , ',kValue(i));
      end
      fprintf(fileID,'R1= %3.5f , R2= %3.5f \n',Rsq1,Rsq2);


      figure(201)
      hold off
      errorbar(LOOHf , Meas , SD , '*k')
      hold on
      errorbar(LOOHf2 , Meas2 , SD2 , '*r')
      plot(x(:,1) , Predict , '-r')
      xlabel('[LOOH]_f (\mu M)')
      ylabel('Crypt survival (%)')
      grid on
      drawnow

end
