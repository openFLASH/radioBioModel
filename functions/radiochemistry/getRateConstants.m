%% getRateConstants
% Return the rate constant for the model
%
%% Syntax
% |[kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants()|
%
%
%% Description
% |[kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants()| Description
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
% |kLOOself| - _SCALAR VECTOR_ - %TERMINATION ROO* + ROO* diffusion limited . Diffusion in lipid membrane
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


function [kbr2 , kLOOself , kb3 , kb8 , kbr] = getRateConstants()


  %TERMINATION R* + R* diffusion limited . Diffusion in lipid membrane
  kbr2 =210162.71223 ;

  %TERMINATION ROO* + ROO* diffusion limited . Diffusion in lipid membrane
  kLOOself =   11123035.47267;

  %PROPAGATION R* + O2
  %kb3 = 5e7  ; %r* + O2 Reference : [1]
  %kb3 = 3e8  ; %r* + O2 Reference : [6]
  kb3 = 182956412.30975  ;

  %TERMINATION LOO* + TOH
  %kb8 = 5.8e3 .* 0.2e-3; %kb8 .* [TOH] Reference :[5]
  %kb8 = 0.04; %Reference : [1]
  kb8 = 0.04522;

  %TERMINATION L* + scav -> LH
  kbr = 119.21460 ;


end
