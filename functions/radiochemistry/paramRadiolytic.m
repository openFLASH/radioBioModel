%% paramRadiolytic
% Define the radiolytic parameters
%
%% Syntax
% |param = paramRadiolytic()|
%
%
%% Description
% |param = paramRadiolytic()| Description
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
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)
%
function param = paramRadiolytic()

%Define the raiolysis rate Constants
%===================================
param.NbSpecies = 9; %Nb of chemical species to track

% --IN WATER---

% These G-value from [7] are for low LET radiation such as 60Co γ-rays and the fast electrons.
ge = 2.8; % nb of aqueous electron per 100 eV [7]
gOHr = 2.8; %nb of radicals OH^. per 100 eV [7]
gHr = 0.62; %nb of radicals H^. per 100 eV [15]
gH2 = 0.47; %nb of H2 per 100 eV [15]
gH2O2 = 0.73; %nb of H202 per 100 eV [7]
%gHp = 2.8; %nb of H+ per 100 eV [7]
% Note that the g-value reported in figure 2 of [14] at 1us post irradiation for proton and electron of various LET are of the same order of magnitude

% --- DIRECT IONISATION IN CELLS -----
gR = 0.59; %Nb of DNA or lipid molecules directly ionised per 100 eV of the radiation

param.kre = ge2kr(ge); %mol/l/Gy Radiolysis rate constant
param.krOHr = ge2kr(gOHr); %mol/l/Gy Radiolysis rate constant
param.krHr = ge2kr(gHr); %mol/l/Gy Radiolysis rate constant
param.krH2 = ge2kr(gH2); %mol/l/Gy Radiolysis rate constant
param.krH2O2 = ge2kr(gH2O2); %mol/l/Gy Radiolysis rate constant
%param.krHp = ge2kr(gHp); %mol/l/Gy Radiolysis rate constant
param.krR = ge2kr(gR); %mol/l/Gy Radiolysis rate constant

%param.kROOself = 1e6;  % [13] rate constant for self reaction of ROO. radical: range from $10^5$ to $10^9 (mol/l)^{-1}s^{-1}$
param.kROOself = 1e4;
param.kde = 1.4e8 .* 1; %k*[DNA] ~ k*[RH] Rate constant for reaction of aqueous electron with bioloigcal molecules (s^-1) [3]
param.kdH = 1e8 ; %TODO Rate constant for reaction of H^. with bioloigcal molecules (s^-1)
param.kdOH = 1e10 .* 6.5e-3; % Rate constant for reaction of HO^. with thiols
% [1] k = 1e10 M-1 s-1 for reaction of OH^. with thiols SH
% [2] GSH concentration 6.5mM



end
