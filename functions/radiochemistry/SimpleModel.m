%-------------------------------------------------
% The simple ODE system
%-------------------------------------------------
function [dydt , labels] = SimpleModel(t,C, k, O2 , DRp , tin)

  labels = {};

  if nargin <1
    dydt=0;
    labels = {'OH^.','L^.','LOO^.' , 'LOOH'};
    return
  end

  k = abs(k);
  DRp_t = interp1(tin,DRp,t,'linear',0);

  %Definition of rate constants
  g_OH = k(1);
  kscv1 = k(2); %OH* scavenger
  kscv2 = k(3); %L* scavenger
  kscv3 = k(8); %LOO* scavenger
  k1 = k(4); % k1 .*  [LH]cst
  k2 = k(5); % L* + L*
  k3 = k(6); % L* + O2
  k4 = k(7); %LOO*

  %Definition of indices
  d_OHr   = 1;
  d_Lr    = 2;
  d_LOOr  = 3;
  d_LOOH  = 4;

  %Concentration at time t
  OHr = C(d_OHr);
  Lr  = C(d_Lr);
  LOOr = C(d_LOOr);
  LOOH = C(d_LOOH);

  %ODE
  dydt = zeros(numel(C),1); %The rate vector
  dydt(d_OHr) = g_OH .* DRp_t - kscv1 .* OHr - k1 .* OHr;
  dydt(d_Lr) = k1 .* OHr - kscv2 .* Lr - k2 .* Lr.^2 - k3 .* Lr .* O2;
  dydt(d_LOOr) = k3 .* Lr .* O2 - k4 .* LOOr - kscv3 .* LOOr;
  dydt(d_LOOH) = k4 .* LOOr;

end
