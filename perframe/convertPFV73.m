function convertPFV73(expdir)

dd = dir(fullfile(expdir,'perframe','*.mat'));

fprintf('Converting %d perframe files\n',numel(dd));  

H = waitbar(0,'converting perframe features..');
for ndx = 1:numel(dd)
  pfname = fullfile(expdir,'perframe',dd(ndx).name);
  ss = fopen(pfname,'r');
  qq = fread(ss,10,'*char')';
  fclose(ss);
  if ~strcmp(qq(8:10),'7.3')
    Q = load(pfname);
    save(pfname,'-struct','Q','-v7.3');
  end
  waitbar(ndx/numel(dd),H);
end
close(H)