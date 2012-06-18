classifierparams = ReadClassifierParamsFile('../classifierparams.txt');
fid = fopen('TrainingDataInfo.txt','w');
for i = 1:numel(classifierparams),
  classifierfile = classifierparams(i).classifierfile;
  classifier = load(classifierfile);
  behavior = classifierparams(i).behaviors.names;
  if iscell(behavior),
    behavior = sprintf('%s_',behavior{:});
    behavior = behavior(1:end-1);
  end
  fprintf('%s:\n',behavior); %#ok<PFCEL>
  nnegbouts_total = 0;
  nposbouts_total = 0;
  nposframes_total = 0;
  nnegframes_total = 0;
  fprintf(fid,'%s\n\n',behavior);
  for j = 1:numel(classifier.expdirs),
    idxneg = strcmpi(classifier.trainingdata(j).names,'None');
    nnegbouts = nnz(idxneg);
    nposbouts = numel(classifier.trainingdata(j).names) - nnegbouts;
    nposframes = sum(classifier.trainingdata(j).t1s(~idxneg) - classifier.trainingdata(j).t0s(~idxneg) + 1);
    nnegframes = sum(classifier.trainingdata(j).t1s(idxneg) - classifier.trainingdata(j).t0s(idxneg) + 1);
    nposbouts_total = nposbouts_total + nposbouts;
    nnegbouts_total = nnegbouts_total + nnegbouts;
    nposframes_total = nposframes_total + nposframes;
    nnegframes_total = nnegframes_total + nnegframes;
    
    fprintf('%s: npositive bouts = %d, nnegative bouts = %d, npos frames = %d, nneg frames = %d\n',classifier.expnames{j},...
      nposbouts,nnegbouts,nposframes,nnegframes);
    fprintf(fid,'%s,%d,%d,%d,%d,%d\n',classifier.expnames{j},classifier.expnames{j}(1)=='p',nposframes,nnegframes,nposbouts,nnegbouts);
  end
  fprintf(fid,'\n');
  
  fprintf('In total for %s, npositive bouts = %d, nnegative bouts = %d, npos frames = %d, nneg frames = %d\n',behavior,...
      nposbouts_total,nnegbouts_total,nposframes_total,nnegframes_total);
  
end
fclose(fid);
% 
% Stop:
% pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110805T101647: npositive bouts = 72, nnegative bouts = 127, npos frames = 760, nneg frames = 1892
% In total for Stop, npositive bouts = 72, nnegative bouts = 127, npos frames = 760, nneg frames = 1892
% WingGrooming:
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 8, nnegative bouts = 25, npos frames = 397, nneg frames = 1526
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 32, nnegative bouts = 74, npos frames = 1278, nneg frames = 3267
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% EXT_CSMH_None_Rig1Plate15BowlA_20120519T170213: npositive bouts = 4, nnegative bouts = 15, npos frames = 234, nneg frames = 1112
% FCF_attP2_1500062_None_Rig1Plate15BowlA_20120519T172815: npositive bouts = 11, nnegative bouts = 21, npos frames = 697, nneg frames = 1390
% FCF_cantons_1500002_None_Rig1Plate15BowlA_20120519T160453: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% pBDPGAL4U_None_Rig1Plate15BowlA_20120519T163429: npositive bouts = 7, nnegative bouts = 15, npos frames = 420, nneg frames = 1298
% In total for WingGrooming, npositive bouts = 62, nnegative bouts = 150, npos frames = 3026, nneg frames = 8593
% Righting:
% pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110610T153218: npositive bouts = 53, nnegative bouts = 310, npos frames = 627, nneg frames = 6093
% pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440: npositive bouts = 21, nnegative bouts = 88, npos frames = 249, nneg frames = 1484
% GMR_14C08_AE_01_TrpA_Rig1Plate15BowlB_20110914T113113: npositive bouts = 27, nnegative bouts = 198, npos frames = 425, nneg frames = 2986
% GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934: npositive bouts = 18, nnegative bouts = 214, npos frames = 146, nneg frames = 4675
% GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039: npositive bouts = 74, nnegative bouts = 180, npos frames = 777, nneg frames = 2407
% In total for Righting, npositive bouts = 193, nnegative bouts = 990, npos frames = 2224, nneg frames = 17645
% Jump:
% GMR_20B07_AE_01_TrpA_Rig2Plate17BowlC_20120216T143940: npositive bouts = 1, nnegative bouts = 57, npos frames = 2, nneg frames = 311
% GMR_14C08_AE_01_TrpA_Rig1Plate15BowlB_20110914T113113: npositive bouts = 19, nnegative bouts = 84, npos frames = 64, nneg frames = 759
% GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934: npositive bouts = 7, nnegative bouts = 78, npos frames = 20, nneg frames = 1038
% GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039: npositive bouts = 193, nnegative bouts = 385, npos frames = 555, nneg frames = 4131
% pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110610T153218: npositive bouts = 0, nnegative bouts = 3, npos frames = 0, nneg frames = 42
% pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440: npositive bouts = 29, nnegative bouts = 61, npos frames = 101, nneg frames = 1385
% FCF_attP2_1500062_None_Rig2Plate17BowlC_20120520T152549: npositive bouts = 20, nnegative bouts = 198, npos frames = 60, nneg frames = 2751
% In total for Jump, npositive bouts = 269, nnegative bouts = 866, npos frames = 802, nneg frames = 10417
% pivot_tail:
% GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039: npositive bouts = 107, nnegative bouts = 132, npos frames = 498, nneg frames = 1166
% GMR_80A01_AE_01_TrpA_Rig2Plate14BowlD_20110408T140618: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 138, nnegative bouts = 155, npos frames = 627, nneg frames = 1306
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 90, nnegative bouts = 227, npos frames = 417, nneg frames = 2294
% FCF_attP2_1500062_None_Rig2Plate17BowlB_20120519T165213: npositive bouts = 19, nnegative bouts = 45, npos frames = 83, nneg frames = 325
% In total for pivot_tail, npositive bouts = 354, nnegative bouts = 559, npos frames = 1625, nneg frames = 5091
% Backup:
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 24, nnegative bouts = 65, npos frames = 165, nneg frames = 2125
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 6, nnegative bouts = 33, npos frames = 30, nneg frames = 967
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804: npositive bouts = 6, nnegative bouts = 3, npos frames = 40, nneg frames = 78
% GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% In total for Backup, npositive bouts = 36, nnegative bouts = 101, npos frames = 235, nneg frames = 3170
% Crabwalk:
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 48, nnegative bouts = 43, npos frames = 464, nneg frames = 1149
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 13, nnegative bouts = 7, npos frames = 204, nneg frames = 104
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804: npositive bouts = 11, nnegative bouts = 20, npos frames = 130, nneg frames = 383
% GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818: npositive bouts = 1, nnegative bouts = 0, npos frames = 4, nneg frames = 0
% In total for Crabwalk, npositive bouts = 73, nnegative bouts = 70, npos frames = 802, nneg frames = 1636
% Walk:
% pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110805T101647: npositive bouts = 156, nnegative bouts = 124, npos frames = 1934, nneg frames = 1138
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804: npositive bouts = 11, nnegative bouts = 19, npos frames = 78, nneg frames = 133
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 28, nnegative bouts = 25, npos frames = 102, nneg frames = 168
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% GMR_80A01_AE_01_TrpA_Rig2Plate14BowlD_20110408T140618: npositive bouts = 16, nnegative bouts = 1, npos frames = 530, nneg frames = 8
% GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039: npositive bouts = 15, nnegative bouts = 54, npos frames = 99, nneg frames = 428
% In total for Walk, npositive bouts = 226, nnegative bouts = 223, npos frames = 2743, nneg frames = 1875
% Chase:
% GMR_14C07_AE_01_TrpA_Rig1Plate15BowlA_20120404T141155: npositive bouts = 78, nnegative bouts = 189, npos frames = 1588, nneg frames = 4308
% GMR_53B02_AE_01_TrpA_Rig1Plate15BowlC_20110930T140347: npositive bouts = 3, nnegative bouts = 8, npos frames = 41, nneg frames = 223
% GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934: npositive bouts = 46, nnegative bouts = 57, npos frames = 821, nneg frames = 1174
% GMR_94B10_AE_01_TrpA_Rig1Plate15BowlC_20111007T155325: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 8, nnegative bouts = 12, npos frames = 124, nneg frames = 264
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 7, nnegative bouts = 8, npos frames = 95, nneg frames = 134
% In total for Chase, npositive bouts = 142, nnegative bouts = 274, npos frames = 2669, nneg frames = 6103
% Touch:
% GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818: npositive bouts = 38, nnegative bouts = 57, npos frames = 336, nneg frames = 1655
% pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928: npositive bouts = 54, nnegative bouts = 77, npos frames = 574, nneg frames = 2735
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545: npositive bouts = 17, nnegative bouts = 25, npos frames = 63, nneg frames = 371
% pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804: npositive bouts = 0, nnegative bouts = 5, npos frames = 0, nneg frames = 45
% EXT_CSMH_None_Rig1Plate15BowlA_20120519T170213: npositive bouts = 0, nnegative bouts = 3, npos frames = 0, nneg frames = 10
% FCF_attP2_1500062_None_Rig1Plate15BowlA_20120519T172815: npositive bouts = 0, nnegative bouts = 2, npos frames = 0, nneg frames = 68
% FCF_cantons_1500002_None_Rig1Plate15BowlA_20120519T160453: npositive bouts = 0, nnegative bouts = 1, npos frames = 0, nneg frames = 4
% pBDPGAL4U_None_Rig1Plate15BowlA_20120519T163429: npositive bouts = 0, nnegative bouts = 0, npos frames = 0, nneg frames = 0
% In total for Touch, npositive bouts = 109, nnegative bouts = 170, npos frames = 973, nneg frames = 4888