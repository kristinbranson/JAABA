function label = readLabels(labelname)
    version = 0.1;

    fid = fopen(labelname,'rb');
    if fid < 0,
      error('Could not open file %s for reading',labelname);
    end

    % format:
    % version (double)
    % moviename_length (int)
    % moviename (char*)
    % matname_length (int)
    % matname (char*)
    % trxname_length (int)
    % trxname (char*)
    % nflies (int)
    % fly ids (int*)
    % firstframe (int)
    % lastframe (int)
    % nbehaviors (int)
    % behaviorname1_length (int)
    % behaviorname1 (char*)
    % ...
    % behaviornamen_length (int)
    % behaviornamen (char*)
    % nbouts (int)
    % startframe_1 (int)
    % endframe_1 (int)
    % behavioridx_1 (int)
    % ...
    % startframe_nbouts (int)
    % endframe_nbouts (int)
    % behavioridx_nbouts (int)

    % header
    v = fread(fid,1,'double'); 
    if v ~= version, 
      error('Version read %f does not match readLabels version %f.',v,version);
    end

    moviename = freadstring(fid);
    matname = freadstring(fid);
    trxname = freadstring(fid);
    nflies = double(fread(fid,1,'int'));
    flies = double(fread(fid,nflies,'int'));
    t0 = double(fread(fid,1,'int'))+1;
    t1 = double(fread(fid,1,'int'))+1;
    nbehaviors = double(fread(fid,1,'int'));

    % behaviors
    behaviors = cell(1,nbehaviors);
    for i = 1:nbehaviors,
      behaviors{i} = freadstring(fid);
    end

    % offset the behavior labels?
    LABELIDXOFFSET = -1;
    if ~ismember('*Unknown*',behaviors),
      LABELIDXOFFSET = LABELIDXOFFSET + 1;
    end
    if ~ismember('*Tracker Failure*',behaviors),
      LABELIDXOFFSET = LABELIDXOFFSET + 1;
    end

    % number of bouts
    nbouts = double(fread(fid,1,'int'));

    % read each bout
    data = double(reshape(fread(fid,nbouts*3,'int'),[3,nbouts]));
    segstarts = data(1,:);
    segends = data(2,:);
    labelidx = data(3,:); % - LABELIDXOFFSET;
    labels = behaviors(labelidx+1);

    fclose(fid);

    label.segstarts = segstarts;
    label.segends = segends;
    label.labels = labels;
    label.behaviors = behaviors;
    label.moviename = moviename;
    label.matname = matname;
    label.trxname = trxname;
    label.t0 = t0;
    label.t1 = t1;
    label.flies = flies;
end

function s = freadstring(fid)
    l = double(fread(fid,1,'int'));
    s = char(fread(fid,l,'char')');
end