%% getRateConstants
% Return the rate constant for the model
%
%% Syntax
% |[kbr2 , kROOself , kb3 , kb8 , kbr] = getRateConstants()|
%
%
%% Description
% |[kbr2 , kROOself , kb3 , kb8 , kbr] = getRateConstants()| Description
%
%
%% Input arguments
% None
%
%
%% Output arguments
%
% |kbr2| - _SCALAR VECTOR_ -  TERMINATION R* + R* diffusion limited . Diffusion in lipid membrane
%
% |kROOself| - _SCALAR VECTOR_ - %TERMINATION ROO* + ROO* diffusion limited . Diffusion in lipid membrane
%
% |kb3| - _SCALAR VECTOR_ -  %PROPAGATION R* + O2
%
% |kb8| - _SCALAR VECTOR_ -  %TERMINATION with TOH
%
% |kbr| - _SCALAR VECTOR_ - %TERMINATION L* + scav -> LH
%
%
%% Contributors
% R.Labarbe
%
%
%% References
% [1] Labarbe, R., Hotoiu, L., Barbier, J. & Favaudon, V. A physicochemical model of reaction kinetics supports peroxyl radical recombination as the main determinant of the FLASH effect. Radiother. Oncol. 1–8 (2020). doi:10.1016/j.radonc.2020.06.001
% [5] Salvador, A. et al. Kinetic Modelling of in Vitro Lipid Peroxidation Experiments - ’ Low Level ’ Validation of a Model of in Vivo Lipid Peroxidation KINETIC MODELLING OF IN VZTRO LIPID. Free Rad. Rex 23, 151–172 (1995)
% [6] Antunes, F., Salvador, A., Marinho, H. S., Alves, R. & Ruy E Pinto. Lipid peroxidation in mitochondrial inner membranes . I . An integrative kinetic model. Free Radic. Biol. Med. 21, 917–943 (1996).


function [kbr2 , kROOself , kb3 , kb8 , kbr , keLOO , koxy , koxyD] = getRateConstants()


  %TERMINATION R* + R* diffusion limited . Diffusion in lipid membrane
  kbr2 =126353.342335   ;

  %TERMINATION ROO* + ROO* diffusion limited . Diffusion in lipid membrane
  kROOself = 10161.179900   ;

  %PROPAGATION R* + O2
  %kb3 = 5e7  ; %r* + O2 Reference : [1]
  %kb3 = 3e8  ; %r* + O2 Reference : [6]
  kb3 = 223503646.894740   ;

  %TERMINATION LOO* + TOH
  %kb8 = 5.8e3 .* 0.2e-3; %kb8 .* [TOH] Reference :[5]
  %kb8 = 0.04; %Reference : [1]
  kb8 = 0.055104 ;

  %TERMINATION L* + scav -> LH
  kbr = 127.012632;

  %Oxilipin rate constants
  keLOO = 10.^(-7.940700);
  koxy  = 10.^-3.717529 ;
  koxyD = 0.019756 ; %Rate constant defined by life time of Oxylipin

end
