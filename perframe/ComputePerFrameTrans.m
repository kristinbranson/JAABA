function [x_trans,IDX,ntrans] = ...
  ComputePerFrameTrans(x,trans_types)

if ~exist('trans_types','var'),
  %trans_types = 'all';
  trans_types = uint8(15);
end

%usealltrans = ischar(trans_types) && strcmpi(trans_types,'all');

x_trans = x;
IDX.orig = 1;
IDX.abs = 0;
IDX.flip = 0;
%if usealltrans || ismember('abs',trans_types),
if bitand(2,trans_types),
  x_trans(end+1,:) = abs(x);
  IDX.abs = size(x_trans,1);
end
%if usealltrans || ismember('flip',trans_types),
if bitand(4,trans_types),
  x_trans(end+1,:) = -x;
  IDX.flip = size(x_trans,1);
end
ntrans = size(x_trans,1);
