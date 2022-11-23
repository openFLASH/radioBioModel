%% acidPartition
% Compute the concentration of the [HA] acid and basic [A] partner of an acid/base couple
% as a function of the [H^+] cocentration.
%
% consider the acid - base equilibrium
% HA <=> A- + H+
% with an equilibrium constant:
% Ka = [A-] [H+]/ [HA]
% pKa = - log10(Ka)
%
% The total concentration : Ct = [A-] + [HA]
%
% Therefore:
% Ct = [HA] . Ka/[Hp] + [HA]  => [HA] = Ct / (1+ Ka/[H+])
% Ct = [A-] + [A-] [H+] / Ka => [A-] = Ct / (1+ [H+] / Ka)
%
%% Syntax
% |[HA , A] = acidPartition(Ct , Hp , pKa)|
%
%
%% Description
% |[HA , A] = acidPartition(Ct , Hp , pKa)| compute the partition into acid and base partner of the couple
%
%
%% Input arguments
% |Ct| - _SCALAR_ - Total concentration (mol/l) of the acid/base couple: Ct = [HA] + [A-]
%
% |Hp| - _SCALAR VECTOR_ - Concentration in H+ (proton). Note that pH = -log(Hp)
%
% |pKa| - _SCALAR_ - pKa = - log10(Ka) where Ka is the equilibrium constant (mol/l) of the acid/base couple
%
%% Output arguments
%
% |HA| - _SCALAR VECTOR_ - |HA(i)| Concentration (mol/l) of the acid partner at Hp(i)
%
% |A| - _SCALAR VECTOR_ - |A(i)| Concentration (mol/l) of the basic partner at Hp(i)
%
% |dHA_dCt| - _SCALAR VECTOR_ - Derivative d[HA]/dCt
%
% |dA_dCt| - _SCALAR VECTOR_ - Derivative d[A]/dCt
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [HA , A , dHA_dCt , dA_dCt ] = acidPartition(Ct , Hp , pKa)

  Ka = power(10,-pKa);
  A = Ct ./ (1+ Hp ./ Ka);
  HA = Ct ./ (1+ Ka ./ Hp);

  %Compute the derivative for the Jacobian
  dA_dCt = 1 ./ (1+ Hp ./ Ka);
  dHA_dCt = 1 ./ (1+ Ka ./ Hp);

end
