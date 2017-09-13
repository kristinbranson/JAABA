function doesmatch = lenientMatch(sin,sout) 

if strcmpi(sin,sout),
  doesmatch = true;
  return;
end

rin = sin;
nins = regexp(rin,'(\d+)','tokenExtents');
if ~isempty(nins),
  rinnew = rin(1:nins{1}(1)-1);
  for k = 1:numel(nins),
    rinnew = [rinnew,'(0*',num2str(str2double(rin(nins{k}(1):nins{k}(2)))),')']; %#ok<AGROW>
    if k < numel(nins),
      rinnew = [rinnew,rin(nins{k}(2)+1:nins{k+1}(1)-1)]; %#ok<AGROW>
    else
      rinnew = [rinnew,rin(nins{k}(2)+1:end)]; %#ok<AGROW>
    end
  end
  rin = rinnew;
end
rin = regexprep(rin,'[-_]','\(\[-_\]\)');

k = regexpi(sout,['^',rin,'$']);
doesmatch = ~isempty(k);
