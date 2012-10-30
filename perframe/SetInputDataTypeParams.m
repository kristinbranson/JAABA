function InputDataTypes = SetInputDataTypeParams()

InputDataTypes = struct;

InputDataTypes.Ctrax = struct;
InputDataTypes.Ctrax.name = 'Ctrax';
InputDataTypes.Ctrax.videorequired = true;
InputDataTypes.Ctrax.hasarena = 'maybe';
InputDataTypes.Ctrax.files = [];

file = struct;
file.name = 'Trx mat file';
file.description = 'Mat file containing trajectories output by Ctrax';
file.code = 'intrxfile';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.Ctrax.files = structappend(InputDataTypes.Ctrax.files,file);

file.name = 'Ann file';
file.code = 'annfile';
file.description = 'Annotation file output by Ctrax';
file.required = true;
file.exts = {'*.ann','*.*'};
InputDataTypes.Ctrax.files = structappend(InputDataTypes.Ctrax.files,file);

InputDataTypes.LarvaeRiveraAlba = struct;
InputDataTypes.LarvaeRiveraAlba.name = 'Rivera-Alba''s Larva Tracker';
InputDataTypes.LarvaeRiveraAlba.videorequired = true;
InputDataTypes.LarvaeRiveraAlba.hasarena = 'must';
InputDataTypes.LarvaeRiveraAlba.files = [];

file = struct;
file.name = 'Trx mat file';
file.code = 'intrxfile';
file.description = 'Mat file containing trajectories output by Marta''s larva tracker';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.LarvaeRiveraAlba.files = structappend(InputDataTypes.LarvaeRiveraAlba.files,file);

InputDataTypes.MouseHouse = struct;
InputDataTypes.MouseHouse.name = 'MouseHouse Tracker';
InputDataTypes.MouseHouse.videorequired = true;
InputDataTypes.MouseHouse.hasarena = 'never';
InputDataTypes.MouseHouse.files = [];

file = struct;
file.name = 'Seq index file';
file.code = 'seqindexfile';
file.description = 'Index for video seq file';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.MouseHouse.files = structappend(InputDataTypes.MouseHouse.files,file);

file.name = 'Trx mat file';
file.code = 'intrxfile';
file.description = 'Mat file containing trajectories output by the MouseHouse tracker';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.MouseHouse.files = structappend(InputDataTypes.MouseHouse.files,file);

InputDataTypes.Qtrax = struct;
InputDataTypes.Qtrax.name = 'Qtrax (CADABRA)';
InputDataTypes.Qtrax.videorequired = true;
InputDataTypes.Qtrax.hasarena = 'must_pxpermm_only';
InputDataTypes.Qtrax.files = [];

file = struct;
file.name = 'Feature file';
file.code = 'featfile';
file.description = 'Feature mat file output by Qtrax';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.Qtrax.files = structappend(InputDataTypes.Qtrax.files,file);

file.name = 'ROI mat file';
file.code = 'roifile';
file.description = 'Mat file ROI information';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.Qtrax.files = structappend(InputDataTypes.Qtrax.files,file);

InputDataTypes.MAGATAnalyzer = struct;
InputDataTypes.MAGATAnalyzer.name = 'MAGATAnalyzer';
InputDataTypes.MAGATAnalyzer.videorequired = true;
InputDataTypes.MAGATAnalyzer.hasarena = 'maybe';
InputDataTypes.MAGATAnalyzer.files = [];

file = struct;
file.name = 'Experiment mat file';
file.code = 'expfile';
file.description = 'Mat file output my the MAGATAnalyzer';
file.required = true;
file.exts = {'*.mat','*.*'};
InputDataTypes.MAGATAnalyzer.files = structappend(InputDataTypes.MAGATAnalyzer.files,file);

InputDataTypes.MWT = struct;
InputDataTypes.MWT.name = 'Multi-Worm Tracker';
InputDataTypes.MWT.videorequired = false;
InputDataTypes.MWT.hasarena = 'maybe';
InputDataTypes.MWT.files = [];

file = struct;
file.name = 'Blobs file';
file.code = 'blobsfile';
file.description = 'Blobs file output my the MWT';
file.required = true;
file.exts = {'*.blobs','*.*'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Outline file';
file.code = 'outlinefile';
file.description = 'Outline file output my the MWT';
file.required = false;
file.exts = {'*.outline','*.*'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Spine file';
file.code = 'spinefile';
file.description = 'Spine file output my the MWT';
file.required = false;
file.exts = {'*.spine','*.*'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Dat file';
file.code = 'datfile';
file.description = 'Dat file output my the MWT containing derived statistics';
file.required = false;
file.exts = {'*.dat','*.*'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

