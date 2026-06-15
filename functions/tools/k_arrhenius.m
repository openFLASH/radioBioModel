%% kD_DMPC
% Rate constant for the diffusion of DMPC (1,2-dimyristoyl-sn- glycero-3-phosphocholine)
% as a function of temperature as measured in [1]
%
%     k = k298 .* f(T)
%
%% Syntax
% |k = kD_DMPC(Tin , k298)|
%
%
%% Description
% |k = kD_DMPC(Tin , k298)| Description
%
%
%% Input arguments
% |Tin| - _SCALAR_ - Temperature (°C)
%
% |k298| -_SCALAR_- Value of the rate constant at 25°C (298K)
%
% |Ea| -_SCALAR_- Activation energy (J/mol)
%
%% Output arguments
%
% |k| - _SCALAR_ - Value of the rate constant at the specifed temperature
%
%
%% REFERENCE
% [1] Bag, N., Hui, D., Yap, X. & Wohland, T. Temperature dependence of diffusion in model and live cell membranes characterized by imaging fl uorescence correlation spectroscopy. Biochim. Biophys. Acta 1838, 802–813 (2014).
%
%% Contributors
% Authors : R. LAbarbe (open.reggui@gmail.com)

function k = k_arrhenius(Tin , k298 , Ea)

  k = k298 .* expA(Tin , Ea) ./ expA(25 , Ea);

end

function k = expA(T , Ea)
  R = 8.31446261815324; % J mol−1 K−1
  T = T + 273; %K
  k = exp(-Ea ./ (R .* T));
end
