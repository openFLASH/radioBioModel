%------------------------------
% sigmoid
%------------------------------
function P = sigmoid(param , X)
  L50 = param(1);
  gam = param(2);
  A = param(3);
  base = param(4);
  P = base + A ./ (1 +  exp(-(X-L50)./gam));
end
