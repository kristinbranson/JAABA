function [data,units] = compute_apt_social(trx,n,fn,apt_info)

parts = strsplit(fn,'_');
view_str = regexp(parts{1},'view(\d*)','tokens');
view = str2double(view_str{1}{1});
fn_type = parts{2};
comp_type = parts{3};
pt1 = str2double(parts{4});
npts = apt_info.params.n_classes;

flies = trx.exp2flies{n};
nflies = numel(flies);


pair = false;
if strcmp(fn_type,'socialpair')
  pt2 = str2double(parts{5});
  pair = true;
else
  pt2 = 1:npts;
  pair = false;
end

[dist,xy,units] = compute_apt_distclosest(trx,n,apt_info,pt1,pt2);

switch comp_type
  case 'dist'
    data = dist;
    units = parseunits('mm');
  case 'ddist'
    data = {};
    for ndx = 1:numel(dist)
      curd = diff(dist{ndx},1,2);
      data{ndx} = curd;
    end
    units = parseunits('mm/s');
  case {'sin','cos'}
    data = {};
    for ndx = 1:numel(dist)      
      xyangle = atan2(xy{ndx}(2,:),xy{ndx}(1,:));
      dangle = xyangle - trx(flies(ndx)).theta;
      if strcmp(comp_type,'sin')
        data{ndx} = sin(dangle);
      else
        data{ndx} = cos(dangle);
      end
    end

end


