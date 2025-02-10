%% radiolysisKinetics
% Reaction rate for the kinetics of the water radiolysis reactions at time |t|.
% This function defines a system of ordinary differnetial equation that can be used as input to the ode15s function.
%
%
%% Syntax
% |dydt = radiolysisKinetics(t,C,param,TissueParam)|
%
%
%% Description
% |dydt = radiolysisKinetics(t,C,param,TissueParam)| Compute radiolysis kinetics
%
%
%% Input arguments
% |t| - _SCALAR_ - time (s)
%
% |C| - _SCALAR MATRIX_ - |C(sp)| Concentrations (u-mol/l) at time |t| for the species |sp|.
% The concentration are specified in u-mol/l so that the tolerance of the ode functions expressed as fractrion of micro-moles
%  *|C(1)| : aqueous electron
%  *|C(2)| : O2 Oxygen
%  *|C(3)| : Ct = [H2O2] + [HO2m]
%  *|C(4)| : Ct = [OHr] + [Orm]
%  *|C(5)| : Hr
%  *|C(6)| : H2
%  *|C(7)| : Ct = [HO2r] + [O2rm] superoxyde radical
%  *|C(8)| : Rr alkyl radicals
%  *|C(9)| : ROO. peroxyl radicals
%
% |param| - _STRUCTURE_ - Radio-chemical kinetic constants
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
% * |param.Fe2p| -_SCALAR_- (uM) [OPTIONAL, default = 0.8873uM] Intracellular concentration of labile Fe++
% * |param.kb8| -_SCALAR_- [OPTIONAL, default = 0.0408 s-1] Rate constant for decay of ROOr
%
% |TissueParam| - _STRUCTURE_ - Data relzated to the tissue
% * |TissueParam.D| -_SCALAR_- m2/s Molecular diffusion constant
% * |TissueParam.pH| -_SCALAR_- pH of the extra vascular tissue. %NB: Assumed buffered solution => ignore production of H+
%
%% Output arguments
%
% |dydt| - _SCALAR MATRIX_ - |dydt(sp,x)|  Reaction rate (u-mol/l/s) of substances |sp| at position |x|
%
% |labels| - _CELL VECTOR_ - |labels{sp}| Name of the chemical species for hte i-th rate |dydt(sp)|
%
% |Jac| -_SCALAR MATRIX_-  Jacobian matrix Jac(sp,j) = d_dydt(sp) / d_C(j) at time t, i.e. derivative of the reaction rate of species sp with respect to species j
%
%% REFERENCES
%
% [1] Azzam, E. I., Jay-Gerin, J.-P., & Pain, D. (2012). Ionizing radiation-induced metabolic oxidative stress and prolonged cell injury. Cancer Letters, 327(1–2), 48–60. http://doi.org/10.1016/j.canlet.2011.12.012
% [2] Buxton, G. V, Greenstock, C. L., Helman, W. P., Ross, A. B. (1988). Critical Review of Rate Constants for Reactions of Hydrated Electrons, Hydrogen Atoms and Hydroxyl Radicals (OH/O−) in Aqueous Solution. Atomic Energy, 17, 513–886. http://doi.org/10.1063/1.555805
% [3] Introduction: radiolysis. (n.d.). http://doi.org/10.1007/s00024-009-0507-0
% [4] Le Caër, S. (2011). Water Radiolysis: Influence of Oxide Surfaces on H2 Production under Ionizing Radiation. Water, 3(4), 235–253. http://doi.org/10.3390/w3010235
% [5] Spitz, D. R., Buettner, G. R., Petronek, M. S., St-aubin, J. J., Flynn, R. T., Waldron, T. J., & Limoli, C. L. (2019). An integrated physico-chemical approach for explaining the differential impact of FLASH versus conventional dose rate irradiation on cancer and normal tissue responses. Radiotherapy and Oncology, (xxxx), 1–5. https://doi.org/10.1016/j.radonc.2019.03.028
% [6] Gray, B., & Carmichael, A. J. (1992). Kinetics of superoxide scavenging by dismutase enzymes and manganese mimics determined by electron spin resonance. Biochemical Journal, 281, 795–802. Retrieved from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1130760/pdf/biochemj00142-0212.pdf
% [7] Neta, P., Robert E. Huie, & Ross, A. B. (1990). Rate Constants for Reactions of Peroxyl Radicals in Fluid Solutions. J. Phys. Chern. Ref. Data, 19, 413–513. https://doi.org/10.1063/1.555854
% [8] https://www.ncbi.nlm.nih.gov/pubmed/22900636
% [9] Reactivity of the hydroxyl radical in aqueous solutions. (1974). National standard reference data system. Retrieved from https://nvlpubs.nist.gov/nistpubs/Legacy/NSRDS/nbsnsrds46.pdf
% [10] QIAN, S. Y., & BUETTNER, G. R. (1999). IRON AND DIOXYGEN CHEMISTRY IS AN IMPORTANT ROUTE TO INITIATION OF BIOLOGICAL FREE RADICAL OXIDATIONS : AN ELECTRON PARAMAGNETIC RESONANCE SPIN TRAPPING STUDY. Free Radical Biology & Medicine, 26, 1447–1456. Retrieved from
% [11] PONKA, P. (1999). Cellular iron metabolism. Kidney International, 55, S2–S11. https://doi.org/10.1046/j.1523-1755.1999.055Suppl.69002.x
% [12] Phaniendra, A., & Babu, D. (2015). Free Radicals : Properties , Sources , Targets , and Their Implication in Various Diseases. Ind J Clin Biochem (Jan-Mar, 30(1), 11–26. https://doi.org/10.1007/s12291-014-0446-0
% [13] Neta, P., Robert E. Huie, & Ross, A. B. (1990). Rate Constants for Reactions of Peroxyl Radicals in Fluid Solutions. J. Phys. Chern. Ref. Data, 19, 413–513. https://doi.org/10.1063/1.555854
% [14] Kesavan, V., Puri, S., & Mohanty, B. (2018). Determination of Kinetics of Peroxidase Enzyme Isolated from Brassica oleracea. IMedPub Journals, 2(2), 2–5. https://doi.org/10.21767/2573-4466.1000
% [15] Britton, B. Y. (1943). The kinetics of the enzyme-substrate of peroxidase. J. Biol. Chem., 151(8), 553–578. Retrieved from http://www.imedpub.com/articles/determination-of-kinetics-of-peroxidase-enzyme-isolated-from-brassica-oleracea.pdf
% [16] Ng, C. F., Schafer, F. Q., Buettner, G. R., & Rodgers, V. G. J. (2007). The rate of cellular hydrogen peroxide removal shows dependency on GSH : Mathematical insight into in vivo H 2 O 2 and GPx concentrations. Free Radical Research, 41(November), 1201?1211. https://doi.org/10.1080/10715760701625075
% [17] Jones, P., & SUGGETT, A. (1968). The Catalase-Hydrogen Peroxide System. Biochem J, 110, 617–620. Retrieved from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1187432/pdf/biochemj00712-0014.pdf
% [18] Barry D. Michael 1986 https://www.ncbi.nlm.nih.gov/pubmed/3741350
% [19] Yasuyuki Ogura ARch. Biochem. Biophys.Catalase Activity at High Concentration of Hydrogen Peroxide 1955
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function [dydt , labels , Jac]= radiolysisKinetics(t,C,param,TissueParam,biologyFLAG)

if nargin < 5
  biologyFLAG = true;
end

%Name of the chemical species tracked by the model
if(biologyFLAG)
  labels = {'e^-_aq','O_2','H_2O_2','OH^.','H^.','H_2','O2^{.-}','R^.','ROO^.'};
else
  labels = {'e^-_aq','O_2','H_2O_2','OH^.','H^.','H_2','O2^{.-}'};
end
if nargin <1
  dydt=0;
  return
end

NbSpecies = length(C); %Number of chemical species involved in radio-chemical reactions

% Define the index label for the tracked species
% This will make it easier to read the code
%==============================================
d_e       = 1; % aqueous electron
d_O2      = 2; % Oxygen
d_CtH2O2  = 3; % H2O2
d_CtOHr   = 4; % radical OH^.
d_Hr      = 5; % Radical H^.
d_H2      = 6; % H2
d_CtH02r  = 7; % Total concentration = [O2rm] superoxyde radical + [H02r]
d_Rr   = 8; %Alkyl radicals
d_ROOr = 9; %peroxylRadical radicals


%Reactant concentrations
%=====================================
wNegatif = find(C<0);
C(wNegatif) = 0; %force positive concentrations
C = C.*1e-6; %convert concentration from u-mol/l into mol/l so that the unit of the constant match the unit of the concentration

%Use more sensible variable name to make code easier to read
e       = C(d_e); % aqueous electron
O2      = C(d_O2); % Oxygen
CtH2O2  = C(d_CtH2O2); % H2O2
CtOHr   = C(d_CtOHr); % radical OH^.
Hr      = C(d_Hr); % Radical H^.
H2      = C(d_H2); % H2
CtH02r  = C(d_CtH02r); % Total concentration = [O2rm] superoxyde radical + [H02r]

if(biologyFLAG)
    Rr   = C(d_Rr); %Alkyl radicals
    ROOr = C(d_ROOr); %peroxylRadical radicals
end

water = 55; %mol/l Constant concentration because excess of water
Hp  = power(10,-TissueParam.pH); %[H+] Buffered solution -> constant concentration
OHm = power(10,-14+TissueParam.pH); %[OH-] Buffered solution -> constant concentration

%Equilibirum
%===========
% OHr + OH- <=> Orm + H2O (equation 14 in [2], eq 15 in [3]) pKa = 11.9
% OHr + H2O <=> Orm + H+
% Ka = [Orm] [Hp] / [OHr]
[OHr , Orm, dOHr_dCt , dOrm_dCt] = acidPartition(CtOHr , Hp , 11.9);

%H2O2 + OHm <=> HO2m + H2O (equation 49 in [2]) pKa = 11.7
%H2O2 + H2O <=> HO2m + Hp
% Ka = [HO2m] [Hp] / [H2O2]
[H2O2 , HO2m, dH2O2_dCt , dHO2m_dCt] = acidPartition(CtH2O2 , Hp , 11.7);

%HO2r <=> O2rm + Hp pKa = 4.9 from [1]
% Ka = [O2rm] [Hp]/ [HO2r]
% The total concentration : Ct = [O2rm] + [HO2r]
[HO2r , O2rm, dHO2r_dCt , dO2rm_dCt] = acidPartition(CtH02r , Hp , 4.9);

RH = get_k(param,'RH',1000e-3); %mol/l [5][10] Assume a constant organic concentration in cells (1 M/l default)

%Declare the individual Jacobians for each reaction rate
%================================
d_R1 = zeros(1,NbSpecies); %d_R1(j) = derivative of the reaction rate R1 with respect to species j
d_R2 = zeros(1,NbSpecies);
d_R3 = zeros(1,NbSpecies);
d_R4 = zeros(1,NbSpecies);
d_R5 = zeros(1,NbSpecies);
d_R6 = zeros(1,NbSpecies);
d_R7 = zeros(1,NbSpecies);
d_R8 = zeros(1,NbSpecies);
d_R9 = zeros(1,NbSpecies);
d_R10 = zeros(1,NbSpecies);
d_R11 = zeros(1,NbSpecies);
d_R12 = zeros(1,NbSpecies);
d_R13 = zeros(1,NbSpecies);
d_R14 = zeros(1,NbSpecies);
d_R15 = zeros(1,NbSpecies);
d_R16 = zeros(1,NbSpecies);
d_R17 = zeros(1,NbSpecies);
d_R18 = zeros(1,NbSpecies);
d_R19 = zeros(1,NbSpecies);
d_R20 = zeros(1,NbSpecies);
d_R21 = zeros(1,NbSpecies);
d_R22 = zeros(1,NbSpecies);
d_R23 = zeros(1,NbSpecies);
d_R24 = zeros(1,NbSpecies);
d_R25 = zeros(1,NbSpecies);
d_R26 = zeros(1,NbSpecies);
d_R27 = zeros(1,NbSpecies);
d_R28 = zeros(1,NbSpecies);
d_R29 = zeros(1,NbSpecies);
d_R30 = zeros(1,NbSpecies);
d_R31 = zeros(1,NbSpecies);
d_R32 = zeros(1,NbSpecies);
d_R33 = zeros(1,NbSpecies);

d_Rb1 = zeros(1,NbSpecies);
d_Rb2 = zeros(1,NbSpecies);
d_Rb3 = zeros(1,NbSpecies);
d_Rb4 = zeros(1,NbSpecies);
d_Rb5 = zeros(1,NbSpecies);
d_Rb6 = zeros(1,NbSpecies);
d_Rb7 = zeros(1,NbSpecies);
d_Rb8 = zeros(1,NbSpecies);
d_Rb9 = zeros(1,NbSpecies);
d_Rb10 = zeros(1,NbSpecies);
d_Rb11 = zeros(1,NbSpecies);

d_Rde = zeros(1,NbSpecies);
d_RdH = zeros(1,NbSpecies);
d_RdOH = zeros(1,NbSpecies);
d_Rbr = zeros(1,NbSpecies);
d_Rbr2 = zeros(1,NbSpecies);

% Reactions from table 2 in [2]
%=====================================
%Reactions with aqueous electron
R1        = 19      .* e .* water;
d_R1(d_e) = 19           .* water; % d_R1(d_e) = the derivative of the reaction rate R1 with respect to species e
% This second line will be used to compute the Jacobian of the system of ODE
% The Jacobian will be used to make ode15s compute faster

R2        = 1.1e10 /2  .* e .* e;
d_R2(d_e) = 1.1e10 /2  .* 2 .* e;

R3        = 2.5e10  .* e .* Hr;
d_R3(d_e) = 2.5e10       .* Hr;
d_R3(d_Hr)= 2.5e10  .* e;

R4           = 3e10    .* e .* OHr;
d_R4(d_e)    = 3e10         .* OHr;
d_R4(d_CtOHr)= 3e10    .* e        .* dOHr_dCt; %d_R4/dCt = (d_R4/dOHr) * (dOHr/dCt)

R5           = 2.2e10  .* e .* Orm;
d_R5(d_e)    = 2.2e10       .* Orm;
d_R5(d_CtOHr)= 2.2e10  .* e         .* dOrm_dCt;

R6        = 2.3e10  .* e .* Hp;
d_R6(d_e) = 2.3e10       .* Hp;

R7            = 1.1e10  .* e .* H2O2;
d_R7(d_e)     = 1.1e10       .* H2O2;
d_R7(d_CtH2O2)= 1.1e10  .* e         .* dHO2m_dCt;

R8            = 3.5e9   .* e .* HO2m;
d_R8(d_e)     = 3.5e9        .* HO2m;
d_R8(d_CtH2O2)= 3.5e9   .* e        .* dHO2m_dCt;

R9        = 1.9e10  .* e .* O2;
d_R9(d_e) = 1.9e10       .* O2;
d_R9(d_O2)= 1.9e10  .* e;

R10            = 1.3e10 .* e .* O2rm;
d_R10(d_e)     = 1.3e10      .* O2rm;
d_R10(d_CtH02r)= 1.3e10 .* e        .* dO2rm_dCt;

%Rde = 0;
Rde       = param.kde .* RH .* e; % Reaction of electron with bioloigical molecules
d_Rde(d_e)= param.kde .* RH;


%REaction with radical H.
R11        = 10  .* Hr .* water;
d_R11(d_Hr)= 10        .* water;

R12        = 1.55e10 /2 .* Hr .* Hr;
d_R12(d_Hr)= 1.55e10 /2 .* 2 .* Hr;

R13           = 7e9     .* Hr .* OHr;
d_R13(d_Hr)   = 7e9           .* OHr;
d_R13(d_CtOHr)= 7e9     .* Hr       .* dOHr_dCt;

R14         = 2.2e7   .* Hr .* OHm;
d_R14(d_Hr) = 2.2e7         .* OHm;

R15            = 9e7     .* Hr .* H2O2;
d_R15(d_Hr)    = 9e7           .* H2O2;
d_R15(d_CtH2O2)= 9e7     .* Hr       .* dH2O2_dCt;

R16 = 2.1e10  .* Hr .* O2;
d_R16(d_Hr) = 2.1e10   .* O2;
d_R16(d_O2) = 2.1e10  .* Hr;

R17            = 1e10    .* Hr .* HO2r;
d_R17(d_Hr)    = 1e10          .* HO2r;
d_R17(d_CtH02r)= 1e10    .* Hr        .* dHO2r_dCt;

%RdH = 0;
RdH        = param.kdH .* RH .* Hr; % Reaction of H^. with biological molecules
d_RdH(d_Hr)= param.kdH .* RH;

%REaction with radical OH.
R18         = 1.1e10 /2  .* OHr .* OHr;
d_R18(d_CtOHr)= 1.1e10 /2  .* 2   .* OHr .* dOHr_dCt;

R20           = 4.2e7    .* OHr .* H2;
d_R20(d_CtOHr)= 4.2e7           .* H2 .* dOHr_dCt;
d_R20(d_H2)   = 4.2e7    .* OHr;

R21           = 1.3e10   .* OHr .* OHm;
d_R21(d_CtOHr)= 1.3e10          .* OHm  .* dOHr_dCt;

R22            = 2.7e7    .* OHr .* H2O2;
d_R22(d_CtOHr) = 2.7e7           .* H2O2 .* dOHr_dCt;
d_R22(d_CtH2O2)= 2.7e7    .* OHr         .* dH2O2_dCt;

R23            = 7.5e9    .* OHr .* HO2m;
d_R23(d_CtOHr) = 7.5e9           .* HO2m .* dOHr_dCt;
d_R23(d_CtH2O2)= 7.5e9    .* OHr         .* dHO2m_dCt;

R25            = 6e9      .* OHr .* HO2r;
d_R25(d_CtOHr) = 6e9             .* HO2r .* dOHr_dCt;
d_R25(d_CtH02r)= 6e9      .* OHr         .* dHO2r_dCt;

R26            = 8e9      .* OHr .* O2rm;
d_R26(d_CtOHr) = 8e9             .* O2rm .* dOHr_dCt;
d_R26(d_CtH02r)= 8e9      .* OHr         .* dO2rm_dCt;

RdOH           = param.kdOH .* OHr; % Reaction of HO^. with thiols
d_RdOH(d_CtOHr)= param.kdOH       .* dOHr_dCt;

%REaction with radical O.-
R27           =  1.8e6  .* Orm .* water;
d_R27(d_CtOHr)=  1.8e6         .* water .* dOrm_dCt;

R32         = 3.6e9   .* Orm .* O2;
d_R32(d_CtOHr)= 3.6e9          .* O2 .* dOrm_dCt;
d_R32(d_O2) = 3.6e9   .* Orm;

R33            = 6e8     .* Orm .* O2rm;
d_R33(d_CtOHr) = 6e8            .* O2rm .* dOrm_dCt;
d_R33(d_CtH02r)= 6e8     .* Orm         .* dO2rm_dCt;

%Reaction with biological molecules
%===================================
if(biologyFLAG)
    %DISMUATION of superoxyde
    % 2H+ + 2 O2rm -> H2O2 + O2
    kb1 = get_k(param,'kb1',2e9);
    kb1p = get_k(param,'kb1p',2e5);
    Rb1            =  kb1 .* O2rm .* O2rm + kb1p .* O2rm .* O2rm; %dismutation of O2rm with superoxyde dismutase [5] 2 O2rm -> O2 + H2O2 2e6 < k 6e9 M-1s-1 [6] and spontaneous dismuation
    d_Rb1(d_CtH02r)= (kb1 .* 2    .* O2rm + kb1p .* 2    .* O2rm) .* dO2rm_dCt;

    %Formation and decay of alkyl radicals Rr
    % OHr + RH -> H2O + Rr
    % Rr + O2 -> ROOr
    % NB: ROOr = peroxyl radical


    %-----------
    %dRH_dt due to OHr
    kb2 = get_k(param,'kb2',1e9);
    Rb2           = kb2 .* OHr .* RH; %dRH_dt  Reactions leading to carbon centered Radicals Rr [5][7][8] [9]table 3 page 22
    d_Rb2(d_CtOHr)= kb2        .* RH .* dOHr_dCt;
    %-----------

    kb3 = get_k(param,'kb3',5e7);
    Rb3        = kb3 .* Rr .* O2; %10^5 to 10^7 in [18]. 10^8 to 10^10 M-1s-1 [5][7].
                          %  5e6 < k < 5e7 M^-1s^-1 Micahel 1981,1986
    d_Rb3(d_O2)= kb3 .* Rr;
    d_Rb3(d_Rr)= kb3       .* O2;

    kbr = get_k(param,'kbr',300);
    Rbr = kbr .* Rr; % This term represents the decay of alkyl radiacl by other routes than oxygen
                    % For example, with thiols: R. + GSH -> RH + RS.
                    % k= 300 Rate constant from Barry D. Michael [18] for untreated
                    % 100 < k[GSH] < 3000 s^-1 Michael 1986, Sonntag 1987
    d_Rbr(d_Rr) = kbr;

    kbr2 = get_k(param,'kbr2',5e7);
    Rbr2         = kbr2 .* Rr .* Rr; % R^. + R^. --> R-R   1e5 < k < 1e9 M^-1s^-1 REference: Favuadon, personal communication
    d_Rbr2(d_Rr) = kbr2 .* 2 .* Rr;


    %Rb4 = 1 .* O2rm .* H2O2; % O2rm + H2O2 -> O2 + OH- + OHr  Haber-Weiss Reaction [10] can be ignored because k is small
    Rb4 =0;

    % REaction with Fe catalyst
    %----------------------------
    if(isfield(param,'Fe2p'))
      [Rb5 , Rb7 , d_Rb5_dH2O2 , dRb7] = ironRadicalsKinetics(H2O2, O2, RH,param.Fe2p);
      d_Rb5(d_CtH2O2) = d_Rb5_dH2O2 .* dH2O2_dCt;
    else
      [Rb5 , Rb7 , d_Rb5_dH2O2 , dRb7] = ironRadicalsKinetics(H2O2, O2, RH);
      d_Rb5(d_CtH2O2) = d_Rb5_dH2O2 .* dH2O2_dCt;
    end

    %Rb5 = dH2O2_dt
    %Rb7 = dRH_dt
    kROOself = get_k(param,'kROOself',1e4);
    Rb6          = kROOself .* ROOr  .* ROOr; %dROOr_dt [13] for self reaction self reaction range from $10^5$ to $10^9 (mol/l)^{-1}s^{-1}$
    d_Rb6(d_ROOr)= kROOself .* 2     .* ROOr;

    kb8 = get_k(param,'kb8',0.0408);
    Rb8          = kb8 .* ROOr;
    d_Rb8(d_ROOr)= kb8 ;

    lipid = get_k(param,'lipid',1e-6); %mol/l
    kb11 = get_k(param,'kb11',20); %(mol/l/s)^(-1)
    Rb11           = kb11 .* lipid .* ROOr; %Chain initiation for unsaturated lipids.
    d_Rb11(d_ROOr) = kb11 .* lipid;

    %Peroxydase reaction
    %-------------------
    % E + S <=Km=> ES --k3--> P
    % AH2 + H2O2 → A + 2 H2O
    %
    Km = 0.44e-6; %mol/l [15] Michaelis constant
    k3 = 3e5; % (mol/l)^-1 . s^-1 [15]
    enzyme = 1e-8; % mol/l Peroxydase concentration in cell [16]
    %Rb9 = k3 .* enzyme.*H2O2 ./ (Km + H2O2); % Michaelis Menten
    Rb9 = 0; %DISABLE PEROXYDASE REACTION

    %Catalase [17]
    %--------------
    % 2 H2O2 → O2 + 2 H2O.
    %https://dial.uclouvain.be/pr/boreal/object/boreal%3A188036/datastream/PDF_01/view
    %
    %Km = 0.9; %mol/l [17] Michaelis constant
    % Km = 1.1; % (M.s)^(-1)
    % k3 = 6.62e7; % (mol/l)^-1 . s^-1 [17]
    % enzyme = 0.08e-6; %mol/l Catalase concentration in cell [19]
    Km = get_k(param,'Km',1.1);% (M.s)^(-1)
    kb10= get_k(param,'kb10',6.62e7); % (mol/l)^-1 . s^-1 [17]
    enzyme = get_k(param,'enzyme',0.08e-6);%mol/l Catalase concentration in cell [19]
    %enzyme = 1e-6.*exp(-t.*0.65); %With the de-activation process [17]
    Rb10           =  kb10.* enzyme.*H2O2 ./ (Km + H2O2);
    d_Rb10(d_CtH2O2) =  dH2O2_dCt .* (kb10 .* enzyme.* (Km + H2O2) - kb10 .* enzyme.*H2O2) ./ (Km + H2O2).^2;

  else
  %Ignore the biological reactions
  Rb1=0;
  Rb2=0;
  Rb3=0;
  Rb4=0;
  Rb5=0;
  Rb6=0;
  Rb7=0;
  Rb8=0;
  Rb9=0;
  Rb10=0;
  Rb11=0;
end %if(biologyFLAG)

%Creation rate of the different species
%=====================================
dydt = zeros(NbSpecies,1); %The rate vector
Jac  = zeros(NbSpecies); %The jacobian


% e  aqueous electron
dydt(d_e) = -R1  - 2.*R2   -R3   -R4   -R5   -R6   -R7   -R8   -R9   -R10   +R14   -Rde + param.kre .* param.R(t,param) ;
Jac(d_e,:)  = -d_R1- 2.*d_R2 -d_R3 -d_R4 -d_R5 -d_R6 -d_R7 -d_R8 -d_R9 -d_R10 +d_R14 -d_Rde;
% The derivative of a sum is equal to the sum of derivative => we are allowed to add the individual Rx
% The Rx are vector. Each element of the vector are added separately. Therefore, each column of the Jacobian matrix gets added separately

% O2 Oxygen
dydt(d_O2) = -R9   -R16   +R26   -R32   +R33   +Rb1   -Rb3   +Rb4   +Rb6   +Rb10;
Jac(d_O2,:)  = -d_R9 -d_R16 +d_R26 -d_R32 +d_R33 +d_Rb1 -d_Rb3 +d_Rb4 +d_Rb6 +d_Rb10;

% Ct = [H2O2] + [HO2m]
dydt(d_CtH2O2) = -R7   -R8   -R15   +R17   +R18   -R23   +Rb1   -Rb4   -Rb5   -Rb9   -2.*Rb10 + param.krH2O2 .* param.R(t,param);
Jac(d_CtH2O2,:)  = -d_R7 -d_R8 -d_R15 +d_R17 +d_R18 -d_R23 +d_Rb1 -d_Rb4 -d_Rb5 -d_Rb9 -2.*d_Rb10;

% Ct = [OHr] + [Orm]
dydt(d_CtOHr) = -R4   -R5   +R7   +R8   +R11   -R13   +R15   -2.*R18   -R20   -R23   -R26   -R32   -R33   -Rb2   +Rb5   -RdOH + param.krOHr .* param.R(t,param);
Jac(d_CtOHr,:)  = -d_R4 -d_R5 +d_R7 +d_R8 +d_R11 -d_R13 +d_R15 -2.*d_R18 -d_R20 -d_R23 -d_R26 -d_R32 -d_R33 -d_Rb2 +d_Rb5 -d_RdOH;
%R21 - R21 : can be ignored: remove OHr and add Orm. Total concentration is unchanged
%-R27 +R27 : can be ignored: remove OHr and add Orm. Total concentration is unchanged

% Radical H^.
dydt(d_Hr) = R1   -R3   +R6   -R11   -2.*R12   -R13   -R14   -R15   -R16   -R17   +R20   -RdH   + param.krHr .* param.R(t,param);
Jac(d_Hr,:)  = d_R1 -d_R3 +d_R6 -d_R11 -2.*d_R12 -d_R13 -d_R14 -d_R15 -d_R16 -d_R17 +d_R20 -d_RdH ;

% H2
dydt(d_H2) = R2   +R3   +R11   +R12   -R20  + param.krH2 .* param.R(t,param);
Jac(d_H2,:)  = d_R2 +d_R3 +d_R11 +d_R12 -d_R20 ;

% Ct = [HO2r] + [O2rm]
dydt(d_CtH02r) = R9   -R10   +R16   -R17   +R23   -R26    -R33   -2.* Rb1   -Rb4 ;
Jac(d_CtH02r,:)  = d_R9 -d_R10 +d_R16 -d_R17 +d_R23 -d_R26  -d_R33 -2.* d_Rb1 -d_Rb4;

if(biologyFLAG)
    dydt(d_Rr) = Rb2   +Rb7   -Rb3   -Rbr   +Rb11   -2 .* Rbr2 + param.krR .* param.R(t,param); %dRr_dt
    Jac(d_Rr,:)  = d_Rb2 +d_Rb7 -d_Rb3 -d_Rbr +d_Rb11 -2 .* d_Rbr2;

    dydt(d_ROOr) = Rb3   -Rb8   -Rb11   -2 .* Rb6; %dROOr_dt [12] for 1st order decay
    Jac(d_ROOr,:)  = d_Rb3 -d_Rb8 -d_Rb11 -2 .* d_Rb6;
end

%Jac % dv/dC = (mol/l.s) / (mol/s) = 1/s => no need to change the concentration units of the Jacobian as the mol/l get cancelled
dydt = dydt .* 1e6; %Convert back  from mol/l into u-mol/l

end
