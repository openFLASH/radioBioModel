%% ge2kr
% Convert the radiation chemical yield Ge (radical / 100eV / incident particle) into the Kr concentation yield (mol/l/Gy)
%
%% Syntax
% |kr = ge2kr(ge)|
%
%
%% Description
% |kr = ge2kr(ge)| convert g into kr
%
%
%% Input arguments
% |ge| - _SCALAR_ - radiation chemical yield Ge (# radical / 100eV / incident particle)
%
%
%% Output arguments
%
% |kr| - _SCALAR_ - Concentration yield (mol/l/Gy)
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function kr = ge2kr(ge)
Na = 6.0221408572e23; %Avogadro's number
q = 1.602176565e-19; %proton charge
d = 1; %kg/l Density of the solvent

C_doseRate = d ./ (q.*100.*Na);
kr = ge .*C_doseRate;

end
