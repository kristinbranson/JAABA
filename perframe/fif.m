function y=fif(test,y_true,y_false)
% "Functional if" --- Like the ternary operator in C.
% Probably only good to use if the 2nd and3rd args are cheap to compute,
% since both will be computed, regardless of the value of test.

if test,
  y=y_true;
else
  y=y_false;
end

end
