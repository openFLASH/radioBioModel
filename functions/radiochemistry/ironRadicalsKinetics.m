%% ironRadicalsKinetics
% Compute the kinetics of the lipid oxidation by H2O2 and O2 catalysed by iron.
%
%% Syntax
% |[Rb5 , Rb7] = ironRadicalsKinetics(H2O2, O2, RH)|
%
%
%% Description
% |[Rb5 , Rb7] = ironRadicalsKinetics(H2O2, O2, RH)| Compute radiolysis kinetics
%
%
%% Input arguments
% |t| - _SCALAR_ - time (s)
%
% |O2| - _SCALAR MATRIX_ - |O2(t)| Concentrations (mol/l) of O2 at time |t|
% |H2O2| -   _SCALAR MATRIX_ - |H2O2(t)| Concentrations (mol/l) of H2O2 at time |t|
% |RH| -   _SCALAR MATRIX_ - |RH(t)| Concentrations (mol/l) of organic material at time |t|
% |Fe2p| -_SCALAR_- (uM) [OPTIONAL, default = 0.8873uM] Intracellular concentration of labile Fe++
%
%
%% Output arguments
%
% |Rb5| - _SCALAR MATRIX_ - |Rb5(t)| Reaction rate of organic material with the H2O2 via Fenton reaction
% |Rb7| - _SCALAR MATRIX_ - |Rb7(t)| Reaction rate of organic material with O2 via Ferryl complexes
%
%% REFERENCES
% [5] Spitz, D. R., Buettner, G. R., Petronek, M. S., St-aubin, J. J., Flynn, R. T., Waldron, T. J., & Limoli, C. L. (2019). An integrated physico-chemical approach for explaining the differential impact of FLASH versus conventional dose rate irradiation on cancer and normal tissue responses. Radiotherapy and Oncology, (xxxx), 1–5. https://doi.org/10.1016/j.radonc.2019.03.028
% [7] Neta, P., Robert E. Huie, & Ross, A. B. (1990). Rate Constants for Reactions of Peroxyl Radicals in Fluid Solutions. J. Phys. Chern. Ref. Data, 19, 413–513. https://doi.org/10.1063/1.555854
% [10] QIAN, S. Y., & BUETTNER, G. R. (1999). IRON AND DIOXYGEN CHEMISTRY IS AN IMPORTANT ROUTE TO INITIATION OF BIOLOGICAL FREE RADICAL OXIDATIONS : AN ELECTRON PARAMAGNETIC RESONANCE SPIN TRAPPING STUDY. Free Radical Biology & Medicine, 26, 1447–1456. Retrieved from
% [10a] Qiang PhD thesis: https://www.healthcare.uiowa.edu/CoreFacilities/esr/education/theses/pdf/QianMSThesisChaps.pdf
% [11] PONKA, P. (1999). Cellular iron metabolism. Kidney International, 55, S2–S11. https://doi.org/10.1046/j.1523-1755.1999.055Suppl.69002.x
% [12] https://pubs.acs.org/doi/10.1021/bi020215g
% [13] Oxygen Free Radicals in Tissue Damage De TARR,M., SAMSON. Springer Science+ Business media LLC. ISBN 978-1-4615-9840-4 https://books.google.be/books?id=XE77AwAAQBAJ&pg=PA34&lpg=PA34&dq=ferryl+complex+with+oxygen&source=bl&ots=izjcNWGv3a&sig=ACfU3U2dWKEWevmvy9AbzCxjHSHBGUiQqQ&hl=fr&sa=X&ved=2ahUKEwiUvYX2j4fjAhUNJVAKHcxrDPAQ6AEwDnoECAcQAQ#v=onepage&q=ferryl%20complex%20with%20oxygen&f=false
%

function [Rb5 , Rb7 , d_Rb5_dH2O2 , dRb7] = ironRadicalsKinetics(H2O2, O2, RH, Fe2p)


% Fenton reaction
%-----------------
if (nargin < 4)
  FeT = 18e-6; %mol/l plasma concentration of Fe [11]
  FeCell =  FeT .* 0.1; % The iron is delivered t the cell via differic transferrin, which represent 10% of the toal amount of Fe-Transferrin [11]
  % In the cell, the iron is in equilibrium betwen ferritin bound and the intracellular iron pool [11]
  Kd = 3 .* 1.48e5; %M^-1 [12] The wild type of the protein is 3 times more efficient than the mutated one in the paper"
  %Kd =  [Fe-Fr] / [Fe2+] .* [FR] ; %dissociation constant of ferritin
  Fe2p = (-1 + sqrt(4.*Kd.*FeCell))./(2.*Kd);
end

kFenton = 1e3; %(Ms)^(-1) %rate constant of the Fenton reaction (about 10^3–10^5 M^-1 s^-1 ) [10]
Rb5        = kFenton .* H2O2 .* Fe2p; %Fe^{2+} + H_2O_2 -> Fe^{3+} + OH^- + OHr Fenton reaction [10]
d_Rb5_dH2O2 = kFenton         .* Fe2p;

% perferryl or ferryl ions oxydations [10]
%-------------------------------------------
%Fe2p + O2 <=> FeO
% Keq = 0.1; %TODO ????
% FeO = Keq .* Fe2p .* O2; %equilibrium constant for Ferryl formation ???? %TODO
% % dRH_dt = kFenton .* Fe2p .* RH
% Rb7 = kFenton .* FeO .* RH; %we assume that the rate constant for oxidation of substrate by ‘Fe2⫹ ⫹ O2’ chemistry is similar to the Fenton reaction ^10$
Rb7=0; % According to [13], these reaction can be ignored.
dRb7 = 0;

end
