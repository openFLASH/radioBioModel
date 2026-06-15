%% find the multiplier of a number

function [c , a , b] = multipliers(y)

  v = 2:y/2;
  c = v(mod(y,v)==0);

  %Find the 2 multipliers close to the square
  if numel(c) > 1
    ysq = sqrt(y);
    [delta, I] = sort(abs(c - ysq));
    a = c(I(1));
    b = c(I(2));
  else
    a = c;
    b = c;
  end

end
