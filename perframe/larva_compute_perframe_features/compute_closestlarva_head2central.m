% closest larva, based on dhead2central
function [data,units,mind] = compute_closestlarva_head2central(trx,n,dosave_d)

 if nargin < 3,
   dosave_d = true;
 end

larvae = trx.exp2flies{n};
nlarvae = numel(larvae);
closestlarva = cell(1,nlarvae);
mind = cell(1,nlarvae);

for i1 = 1:nlarvae,
  larva1 = larvae(i1);
  fprintf('fly1 = %d\n',larva1);
  larvae2 = larvae;
  d = nan(numel(larvae2),trx(larva1).nframes);  
  
  for i2 = 1:numel(larvae2),
    larva2 = larvae2(i2);
    if larva1 == larva2,
      continue;
    end
    d(i2,:) = dhead2central_pair(trx,larva1,larva2);
  end
  [mind{i1},closesti] = min(d,[],1);
  closestlarva{i1} = larvae2(closesti);
  closestlarva{i1}(isnan(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
 if dosave_d,
   data = mind; %#ok<NASGU>
   units = parseunits('mm'); %#ok<NASGU>
   filename = trx.GetPerFrameFile('dhead2central',n);
    
   save(filename,'data','units');
 end

data = closestlarva;
units = parseunits('unit');