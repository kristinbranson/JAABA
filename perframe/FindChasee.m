function FindChasee(explist,chaseJabFile,perframe2use,outjabfilename,blenthres)

Q = load(chaseJabFile,'-mat');
P = load(outjabfilename,'-mat');
outfilename = P.x.file.scorefilename;
for ndx = 1:numel(explist)
  expname = explist{ndx};
  dat = load( fullfile(expname,Q.x.file.perframedir,[perframe2use '.mat']));
  scorefilename = Q.x.file.scorefilename;
  if ~strcmpi(scorefilename(end-3:end),'.mat')
    scorefilename = [scorefilename '.mat'];
  end
  S = load( fullfile(expname,scorefilename));
  
  O = S;
  O.jabFileNameAbs= outjabfilename;
  
  for fly = 1:numel(S.allScores.scores)
      O.allScores.scores{fly} = -ones(size(O.allScores.scores{fly}));
    O.allScores.postprocessed{fly} = zeros(size(O.allScores.postprocessed{fly}));
    O.allScores.t0s{fly} = [];
    O.allScores.t1s{fly} = [];
    O.allScores.scoreNorm = 1;
  end
  
  for fly = 1:numel(S.allScores.scores)
    [t0s t1s] = get_interval_ends(S.allScores.postprocessed{fly}>0);
    tstart = S.allScores.tStart(fly);
    for bndx = 1:numel(t0s)
      curt = t0s(bndx):t1s(bndx)-1;
      chasee = dat.data{fly}( curt -tstart + 1);
      if numel(chasee) ==1
        astart = 1;
        aend = 1;
        blens = 1;
      else
        ss = chasee(2:end)==chasee(1:end-1);
        [astart,aend] = get_interval_ends(ss);
        blens = aend-astart+1;
      end
%       [hh,~,occ] = unique(chasee);
%       counts = histc(occ,1:numel(hh));
%       
%       [a,b] = max(counts);
%       if a/numel(occ)>0.8
%         curidx = hh(b);
%         O.allScores.scores{curidx}(curt) = 1;
%         O.allScores.postprocessed{curidx}(curt) = 1;
%         O.allScores.t0s{curidx}(end+1) = t0s(bndx);
%         O.allScores.t1s{curidx}(end+1) = t1s(bndx);
%       end
      for count = 1:numel(blens)
        if blens(count)>=blenthres
          curidx = chasee(astart(count));
          t = (astart(count):aend(count))+curt(1)-1;
          O.allScores.scores{curidx}(t) = 1;
          O.allScores.postprocessed{curidx}(t) = 1;
          O.allScores.t0s{curidx}(end+1) = astart(count)+curt(1)-1;
          O.allScores.t1s{curidx}(end+1) = aend(count)+curt(1)-1;
        end
      end
    end
    
    
  end
  save(fullfile(expname,outfilename),'-struct','O');
  
end