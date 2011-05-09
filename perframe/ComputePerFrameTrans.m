function [x_trans,IDX_ORIG,IDX_ABS,IDX_FLIP,ntrans] = ...
  ComputePerFrameTrans(x,trans_types)

if ~exist('trans_types','var'),
  trans_types = 'all';
end

usealltrans = ischar(trans_types) && strcmpi(trans_types,'all');

x_trans = x;
IDX_ORIG = 1;
IDX_ABS = 0;
IDX_FLIP = 0;
if usealltrans || ismember('abs',trans_types),
  x_trans(end+1,:) = abs(x);
  IDX_ABS = size(x_trans,1);
end
if usealltrans || ismember('flip',trans_types),
  x_trans(end+1,:) = -x;
  IDX_FLIP = size(x_trans,1);
end
ntrans = size(x_trans,1);