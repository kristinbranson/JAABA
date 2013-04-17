function InputDataTypes = SetInputDataTypeParams()

InputDataTypes = struct;

InputDataTypes.Ctrax = struct;
InputDataTypes.Ctrax.name = 'Ctrax';
InputDataTypes.Ctrax.videorequired = true;
InputDataTypes.Ctrax.readarena = 'maybe';
InputDataTypes.Ctrax.readpxpermm = 'maybe';
InputDataTypes.Ctrax.readfps = 'maybe';
InputDataTypes.Ctrax.writearena = 'maybe';
InputDataTypes.Ctrax.writepxpermm = 'maybe';
InputDataTypes.Ctrax.writefps = 'maybe';
InputDataTypes.Ctrax.files = [];

file = struct;
file.name = 'Trx mat file';
file.description = 'Mat file containing trajectories output by Ctrax';
file.code = 'intrxfile';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.Ctrax.files = structappend(InputDataTypes.Ctrax.files,file);

file.name = 'Ann file';
file.code = 'annfile';
file.description = 'Annotation file output by Ctrax';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.ann'};
InputDataTypes.Ctrax.files = structappend(InputDataTypes.Ctrax.files,file);

InputDataTypes.LarvaeRiveraAlba = struct;
InputDataTypes.LarvaeRiveraAlba.name = 'Rivera-Alba''s Larva Tracker';
InputDataTypes.LarvaeRiveraAlba.videorequired = true;
InputDataTypes.LarvaeRiveraAlba.readarena = 'yes';
InputDataTypes.LarvaeRiveraAlba.readpxpermm = 'yes';
InputDataTypes.LarvaeRiveraAlba.readfps = 'yes';
InputDataTypes.LarvaeRiveraAlba.writearena = 'no';
InputDataTypes.LarvaeRiveraAlba.writepxpermm = 'no';
InputDataTypes.LarvaeRiveraAlba.writefps = 'no';
InputDataTypes.LarvaeRiveraAlba.files = [];

file = struct;
file.name = 'Trx mat file';
file.code = 'intrxfile';
file.description = 'Mat file containing trajectories output by Marta''s larva tracker';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeRiveraAlba.files = structappend(InputDataTypes.LarvaeRiveraAlba.files,file);

InputDataTypes.MoTr = struct;
InputDataTypes.MoTr.name = 'MoTr';
InputDataTypes.MoTr.videorequired = true;
InputDataTypes.MoTr.readarena = 'no';
InputDataTypes.MoTr.readpxpermm = 'no';
InputDataTypes.MoTr.readfps = 'yes';
InputDataTypes.MoTr.writearena = 'yes';
InputDataTypes.MoTr.writepxpermm = 'yes';
InputDataTypes.MoTr.writefps = 'no';
InputDataTypes.MoTr.files = [];

file = struct;
file.name = 'Seq index file';
file.code = 'seqindexfile';
file.description = 'Index for video seq file';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.MoTr.files = structappend(InputDataTypes.MoTr.files,file);

file.name = 'Trx mat file';
file.code = 'intrxfile';
file.description = 'Mat file containing trajectories output by MoTr';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.MoTr.files = structappend(InputDataTypes.MoTr.files,file);

InputDataTypes.Qtrax = struct;
InputDataTypes.Qtrax.name = 'Qtrax (CADABRA)';
InputDataTypes.Qtrax.videorequired = true;
InputDataTypes.Qtrax.readarena = 'no';
InputDataTypes.Qtrax.readpxpermm = 'yes';
InputDataTypes.Qtrax.readfps = 'yes';
InputDataTypes.Qtrax.writearena = 'yes';
InputDataTypes.Qtrax.writepxpermm = 'no';
InputDataTypes.Qtrax.writefps = 'no';
InputDataTypes.Qtrax.files = [];

file = struct;
file.name = 'Feature file';
file.code = 'featfile';
file.description = 'Feature mat file output by Qtrax';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.Qtrax.files = structappend(InputDataTypes.Qtrax.files,file);

file.name = 'ROI mat file';
file.code = 'roifile';
file.description = 'Mat file ROI information';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.Qtrax.files = structappend(InputDataTypes.Qtrax.files,file);

InputDataTypes.MAGATAnalyzer = struct;
InputDataTypes.MAGATAnalyzer.name = 'MAGATAnalyzer';
InputDataTypes.MAGATAnalyzer.videorequired = true;
InputDataTypes.MAGATAnalyzer.readarena = 'no';
InputDataTypes.MAGATAnalyzer.readpxpermm = 'yes';
InputDataTypes.MAGATAnalyzer.readfps = 'yes';
InputDataTypes.MAGATAnalyzer.writearena = 'yes';
InputDataTypes.MAGATAnalyzer.writepxpermm = 'no';
InputDataTypes.MAGATAnalyzer.writefps = 'no';
InputDataTypes.MAGATAnalyzer.files = [];

file = struct;
file.name = 'Experiment mat file';
file.code = 'expfile';
file.description = 'Mat file output by the MAGATAnalyzer';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.MAGATAnalyzer.files = structappend(InputDataTypes.MAGATAnalyzer.files,file);

InputDataTypes.MWT = struct;
InputDataTypes.MWT.name = 'Multi-Worm Tracker';
InputDataTypes.MWT.videorequired = false;
InputDataTypes.MWT.readarena = 'no';
InputDataTypes.MWT.readpxpermm = 'maybe';
InputDataTypes.MWT.readfps = 'yes';
InputDataTypes.MWT.writearena = 'yes';
InputDataTypes.MWT.writepxpermm = 'maybe';
InputDataTypes.MWT.writefps = 'no';
InputDataTypes.MWT.files = [];

file = struct;
file.name = 'Blobs file';
file.code = 'blobsfile';
file.description = 'Blobs file(s) output by the MWT';
file.required = true;
file.multiplefiles = 1;
file.isdir = 0;
file.exts = {'*.blobs'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Spine file';
file.code = 'spinefile';
file.description = 'Spine file output by the MWT';
file.required = false;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.spine'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

file = struct;
file.name = 'Dat files';
file.code = 'datfiles';
file.description = 'Dat file(s) output by the MWT containing derived statistics';
file.required = false;
file.multiplefiles = 1;
file.isdir = 0;
file.exts = {'*.dat'};
InputDataTypes.MWT.files = structappend(InputDataTypes.MWT.files,file);

%% Sam Reid's larvae

InputDataTypes.LarvaeReid = struct;
InputDataTypes.LarvaeReid.name = 'Reid''s Larva Data';
InputDataTypes.LarvaeReid.videorequired = false;
InputDataTypes.LarvaeReid.readarena = 'no';
InputDataTypes.LarvaeReid.readpxpermm = 'maybe';
InputDataTypes.LarvaeReid.readfps = 'yes';
InputDataTypes.LarvaeReid.writearena = 'yes';
InputDataTypes.LarvaeReid.writepxpermm = 'maybe';
InputDataTypes.LarvaeReid.writefps = 'no';
InputDataTypes.LarvaeReid.files = [];

file = struct;
file.name = 'Blobs file';
file.code = 'blobsfile';
file.description = 'Blobs file(s) output by the MWT';
file.required = true;
file.multiplefiles = 1;
file.isdir = 0;
file.exts = {'*.blobs'};
InputDataTypes.LarvaeReid.files = structappend(InputDataTypes.LarvaeReid.files,file);

file = struct;
file.name = 'kinData mat file';
file.code = 'kinmatfile';
file.description = 'kinData mat file computed by Reid''s analysis on top of the MWT output';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeReid.files = structappend(InputDataTypes.LarvaeReid.files,file);

file = struct;
file.name = 'eventData mat file';
file.code = 'eventmatfile';
file.description = 'eventData mat file computed by Reid''s analysis on top of the MWT output';
file.required = false;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeReid.files = structappend(InputDataTypes.LarvaeReid.files,file);

%% Ibraheem's single larva tracker data

InputDataTypes.LarvaeIbraheem = struct;
InputDataTypes.LarvaeIbraheem.name = 'Ibraheem''s Single Larva Tracker Data';
InputDataTypes.LarvaeIbraheem.videorequired = false;
InputDataTypes.LarvaeIbraheem.readarena = 'no';
InputDataTypes.LarvaeIbraheem.readpxpermm = 'maybe';
InputDataTypes.LarvaeIbraheem.readfps = 'no';
InputDataTypes.LarvaeIbraheem.writearena = 'yes';
InputDataTypes.LarvaeIbraheem.writepxpermm = 'maybe';
InputDataTypes.LarvaeIbraheem.writefps = 'no';
InputDataTypes.LarvaeIbraheem.files = [];

file = struct;
file.name = 'First frame image file';
file.code = 'firstframeimfile';
file.description = 'First frame of the video';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.BMP'};
InputDataTypes.LarvaeIbraheem.files = structappend(InputDataTypes.LarvaeIbraheem.files,file);

file = struct;
file.name = 'Cropped images from the video';
file.code = 'croppedimsmatfile';
file.description = 'Mat file containing cropped images around the larva.';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeIbraheem.files = structappend(InputDataTypes.LarvaeIbraheem.files,file);

file = struct;
file.name = 'Bounding box mat file';
file.code = 'bboxmatfile';
file.description = 'Mat file containing bounding boxes of cropped images';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeIbraheem.files = structappend(InputDataTypes.LarvaeIbraheem.files,file);


file = struct;
file.name = 'Contour mat file';
file.code = 'contourmatfile';
file.description = 'Mat file containing larva contours';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeIbraheem.files = structappend(InputDataTypes.LarvaeIbraheem.files,file);


file = struct;
file.name = 'Coordinates mat file';
file.code = 'coordsmatfile';
file.description = 'Mat file containing coordinates of points on the larva, in pixels';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeIbraheem.files = structappend(InputDataTypes.LarvaeIbraheem.files,file);


file = struct;
file.name = 'MotorData mat file';
file.code = 'motormatfile';
file.description = 'Mat file containing motorData';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.LarvaeIbraheem.files = structappend(InputDataTypes.LarvaeIbraheem.files,file);


%%

InputDataTypes.CtraxPlusWings = struct;
InputDataTypes.CtraxPlusWings.name = 'Ctrax plus wing tracker';
InputDataTypes.CtraxPlusWings.videorequired = true;
InputDataTypes.CtraxPlusWings.readarena = 'maybe';
InputDataTypes.CtraxPlusWings.readpxpermm = 'maybe';
InputDataTypes.CtraxPlusWings.readfps = 'maybe';
InputDataTypes.CtraxPlusWings.writearena = 'maybe';
InputDataTypes.CtraxPlusWings.writepxpermm = 'maybe';
InputDataTypes.CtraxPlusWings.writefps = 'maybe';
InputDataTypes.CtraxPlusWings.files = [];

file = struct;
file.name = 'Trx mat file';
file.description = 'Mat file containing trajectories output by Ctrax plus wing tracker';
file.code = 'intrxfile';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.CtraxPlusWings.files = structappend(InputDataTypes.CtraxPlusWings.files,file);

file.name = 'Ann file';
file.code = 'annfile';
file.description = 'Annotation file output by Ctrax';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.ann'};
InputDataTypes.CtraxPlusWings.files = structappend(InputDataTypes.CtraxPlusWings.files,file);

file.name = 'Per-frame directory';
file.code = 'inperframedir';
file.description = 'Per-frame directory containing wing statistics';
file.required = true;
file.multiplefiles = 0;
file.isdir = 1;
file.exts = {'*'};
InputDataTypes.CtraxPlusWings.files = structappend(InputDataTypes.CtraxPlusWings.files,file);

InputDataTypes.SimpleTwoFlies = struct;
InputDataTypes.SimpleTwoFlies.name = 'Kristin''s Simple Two-Fly Tracker';
InputDataTypes.SimpleTwoFlies.videorequired = true;
InputDataTypes.SimpleTwoFlies.readarena = 'no';
InputDataTypes.SimpleTwoFlies.readpxpermm = 'no';
InputDataTypes.SimpleTwoFlies.readfps = 'no';
InputDataTypes.SimpleTwoFlies.writearena = 'maybe';
InputDataTypes.SimpleTwoFlies.writepxpermm = 'yes';
InputDataTypes.SimpleTwoFlies.writefps = 'yes';
InputDataTypes.SimpleTwoFlies.files = [];

file = struct;
file.name = 'Trx mat file';
file.description = 'Mat file containing trajectories and other information output by Kristin''s Simple Two-Fly Tracker';
file.code = 'intrxfile';
file.required = true;
file.multiplefiles = 0;
file.isdir = 0;
file.exts = {'*.mat'};
InputDataTypes.SimpleTwoFlies.files = structappend(InputDataTypes.SimpleTwoFlies.files,file);

file.name = 'Per-frame directory';
file.code = 'inperframedir';
file.description = 'Per-frame directory containing wing statistics';
file.required = true;
file.multiplefiles = 0;
file.isdir = 1;
file.exts = {'*'};
InputDataTypes.SimpleTwoFlies.files = structappend(InputDataTypes.SimpleTwoFlies.files,file);