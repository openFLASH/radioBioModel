%% paramKinectisDiffusion
% Description
%
%% Syntax
% |[param , TissueParam ]= paramKinectisDiffusion()|
%
%
%% Description
% |[param , TissueParam ]= paramKinectisDiffusion()| Description
%
%
%% Input arguments
% None
%
%
%% Output arguments
%
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants
% * |param.NbSpecies| -_INTEGER_- Number of chemical species involved in radio-chemical reactions
% * |param.kre| -_SCALAR_- mol/l/Gy Radiolysis rate constant for generation of aquaous electrons
% * |param.krOHr| -_SCALAR_- mol/l/Gy Radiolysis rate constant for generation of radical OH.
% * |param.krHr| -_SCALAR_- mol/l/Gy Radiolysis rate constant for generation of radical H.
% * |param.krH2| -_SCALAR_- mol/l/Gy Radiolysis rate constant for generation of H2
% * |param.krH2O2| -_SCALAR_- mol/l/Gy Radiolysis rate constant for generation of H2O2
% * |param.krHp| -_SCALAR_- mol/l/Gy Radiolysis rate constant  for generation of H+ NB: Assumed buffered solution => ignore production of H+
% * |param.R| -_FUNCTION POINTER_- Pointer to the function Beam(t) defining if the beam is ON at time t
% * |param.td| -_SCALAR_- (s) Total duration of the sequence of BEAM ON pulses
% * |param.R0| -_SCALAR_- (Gy/s) Dose rate
% * |param.t_on| -_SCALAR_- (s) For pulsed beam, duration of a single pulse
% * |param.T| -_SCALAR_- (s) Period separating 2 pulses
%
% |TissueParam| - _STRUCTURE_ - Data relzated to the tissue
% * |TissueParam.D| -_SCALAR_- m2/s Molecular diffusion constant
% * |TissueParam.pH| -_SCALAR_- pH of the extra vascular tissue. %NB: Assumed buffered solution => ignore production of H+
% * |TissueParam.O2| -_SCALAR_- (uM) Inital oxygen concentration
%
% References
% [1] Held, K. D., Harrop, H. A., & Michel, B. D. (1981). Effects of oxygen and sulphydryl-containing compounds on irradiated transforming DNA . Part I . Actions of dithiothreitol. Int. J. Radiat. Biol, 40(6), 613–622. https://doi.org/10.1080/09553008114551601
% [2] Simic, M. G., Grossman, L., & Upton, A. C. (1986). Mechanisms of DNA damage and repair. In M. G. Simic, L. Grossman, & A. C. Upton (Eds.). London: Plenum Press.
% [3] Sonntag The chemistry of free radical mediated DNA damage Physiocal and chemical mechanism in molecular radiation biology Ed. WA Glass & MN Varma Plenum Press NY 1991
% [4] Cobut, V., Frongillo, Y., Patau, J. P., Goulet, T., Fraser, M. J., Jay-Gerin, J. P. (1998). Monte Carlo simulation of fast electron and proton tracks in liquid water - I. Physical and physicochemical aspects. Radiation Physics and Chemistry, 51(3), 229–243
% [5] Frongillo, Goulet, Fraser, Cobut, Patau and Jay-Gerin Monte carlo simulation of fast electron and proton tracks in liquid water ii. nonhomogeneous chemistry radiat. phys. chem. vol. 51, no. 3, pp. 245-254, 1998
% [7] Kawade, Vitthal Ajinath (n.d.). Synthesis characterization and pulse radiolysis study of mixed ligand cobalt complexes of pyridine and polypyridyl derivatives. 2010. University of Pune. http://hdl.handle.net/10603/2689
% [13] Neta, P., Robert E. Huie, & Ross, A. B. (1990). Rate Constants for Reactions of Peroxyl Radicals in Fluid Solutions. J. Phys. Chern. Ref. Data, 19, 413–513. https://doi.org/10.1063/1.555854
% [14] Frongillo, Y., Goulet, T., Fraser, M.-J., Cobut, V., Patau, J. P., & Jay-Gerin, J.-P. (1998). MONTE CARLO SIMULATION OF FAST ELECTRON AND PROTON TRACKS IN LIQUID WATER—II. NONHOMOGENEOUS CHEMISTRY. Radiation Physics and Chemistry, 51(3), 245–254. https://doi.org/10.1016/S0969-806X(97)00097-2
% [15] 12.	Buxton GV.  An overview of the radiation chemistry of liquids. In: Radiation chemistry, from basics to applications in material and life sciences. Spotheim-Maurizot M, Mostafavi M, Douki T, Belloni J, editors. EDP Sciences, Les Ulis, 2008. pp. 3-16.
% [16] Liset de la Fuente Rosales, A Monte Carlo Study of the Direct and Indirect DNA Damage Induced by Ionizing Radiation, Ph.D. thesis, Universidade Estadual de Campinas (2018).

%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)
%
function [param , TissueParam ]= paramKinectisDiffusion()

%Define the raiolysis rate Constants
%===================================
param = paramRadiolytic();

%Define the diffusion Constants
%==============================
TissueParam.D = 1.7e-5 * 1e-4; %cm^2 /s [7] => m2/s
TissueParam.pH = 7; %pH of the extra vascular tissue.
TissueParam.O2 = 50;% uM Inital oxygen concentration

end
