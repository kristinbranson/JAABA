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
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.Ctrax.files = structappend(InputDataTypes.Ctrax.files,file);

file.name = 'Ann file';
file.code = 'annfile';
file.description = 'Annotation file output by Ctrax';
file.required = true;
file.multiplefiles = 0;
file.exts = {'*.ann'};
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
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeRiveraAlba.files = structappend(InputDataTypes.LarvaeRiveraAlba.files,file);

InputDataTypes.MoTr = struct;
InputDataTypes.MoTr.name = 'MoTr';
InputDataTypes.MoTr.videorequired = true;
InputDataTypes.MoTr.hasarena = 'never';
InputDataTypes.MoTr.files = [];

file = struct;
file.name = 'Seq index file';
file.code = 'seqindexfile';
file.description = 'Index for video seq file';
file.required = true;
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.MoTr.files = structappend(InputDataTypes.MoTr.files,file);

file.name = 'Trx mat file';
file.code = 'intrxfile';
file.description = 'Mat file containing trajectories output by MoTr';
file.required = true;
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.MoTr.files = structappend(InputDataTypes.MoTr.files,file);

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
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.Qtrax.files = structappend(InputDataTypes.Qtrax.files,file);

file.name = 'ROI mat file';
file.code = 'roifile';
file.description = 'Mat file ROI information';
file.required = true;
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.Qtrax.files = structappend(InputDataTypes.Qtrax.files,file);

InputDataTypes.MAGATAnalyzer = struct;
InputDataTypes.MAGATAnalyzer.name = 'MAGATAnalyzer';
InputDataTypes.MAGATAnalyzer.videorequired = true;
InputDataTypes.MAGATAnalyzer.hasarena = 'must_pxpermm_only';
InputDataTypes.MAGATAnalyzer.files = [];

file = struct;
file.name = 'Experiment mat file';
file.code = 'expfile';
file.description = 'Mat file output my the MAGATAnalyzer';
file.required = true;
file.multiplefiles = 0;
file.exts = {'*.mat'};
InputDataTypes.MAGATAnalyzer.files = structappend(InputDataTypes.MAGATAnalyzer.files,file);

InputDataTypes.MWT = struct;
InputDataTypes.MWT.name = 'Multi-Worm Tracker';
InputDataTypes.MWT.videorequired = false;
InputDataTypes.MWT.hasarena = 'must_pxpermm_only';
InputDataTypes.MWT.files = [];

file = struct;
file.name = 'Blobs file';
file.code = 'blobsfile';
file.description = 'Blobs file(s) output my the MWT';
file.required = true;
file.multiplefiles = 1;
file.exts = {'*.blobs'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Spine file';
file.code = 'spinefile';
file.description = 'Spine file output my the MWT';
file.required = false;
file.multiplefiles = 0;
file.exts = {'*.spine'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Dat files';
file.code = 'datfiles';
file.description = 'Dat file(s) output my the MWT containing derived statistics';
file.required = false;
file.multiplefiles = 1;
file.exts = {'*.dat'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

