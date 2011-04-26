%-----------------------------------------------------------------------------
% efdwmaker
%-----------------------------------------------------------------------------
% ExpFun Project
%-----------------------------------------------------------------------------
% Used to pre-process files for subsequent use in ExpFunBrowser program
% 
% Output variables and commands are written to a .efd file 
% (that is in Matlab MAT format!)
%-----------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------
% See Also: expfunbrowser
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Code History
%-----------------------------------------------------------------------------
% Revisions:
%	21 Feb 2011 (SJS):
%	 -	added comments
%	 -	added checks on files
%	19 Apr 2011 (SJS) created efdwmaker.m to deal with datawave data
%-----------------------------------------------------------------------------

%-------------------------------------------------------------
% Get file information and store in files variable
% 	files is a struct array with fields:
% 		name		file name
% 		date		file date information (from system)
% 		bytes		size of file in bytes
% 		isdir		0 if regular file, 1 if directory
% 		datenum	numeric date information
%-------------------------------------------------------------
% first, use dir command to find .src files
files = dir('G1_3_890_RSP1_1175.src');

%-------------------------------------------------------------
% create some struct templates
%-------------------------------------------------------------
% create src struct template
src = struct('noof_files',1,'filename',{''},'date','','size','');
% create plx struct template
plx = struct('noof_files',1,'filename',{''},'startinterval',1,'noof_intervals',-1);

%-------------------------------------------------------------
% store commands to be run in cmds cell array
%-------------------------------------------------------------
% load src file data
cmds{1} = 'srcraster = loadSRCData(srcfiledetails.filename{1});';
% schedule ????
cmds{2} = 'schedule';
% load plx data
cmds{3} = 'plx = loadPLXFile(plx_srcfiledetails.filename{1});';
% convert plx data into rasters
cmds{4} = 'raster = plxIntoRaster(plx,srcraster);';
% generate 5 ms bin raster
cmds{5} = 'psth_5msBW = raster2PSTH(raster,''start'',0,''end'',1.4,''bin_width'',5e-3);';
% ?
cmds{6} = 'psth_5msBW_3500ms = raster2PSTH(raster,''start'',0,''end'',3.5,''bin_width'',5e-3);';
% 1 ms bin raster
cmds{7} = 'psth_1msBW = raster2PSTH(raster,''start'',0,''end'',1.4,''bin_width'',1e-3);';
% generate global 1 ms bin PSTH
cmds{8} = 'globalpsth1msBW_plx = plx2PSTH(plx,''same'',''start'',0,''end'',0.4,''bin_width'',1e-3,''event_chan'',1);';
% generate global 5 ms bin PSTH
cmds{9} = 'globalpsth5msBW_plx = plx2PSTH(plx,''same'',''start'',0,''end'',0.4,''bin_width'',5e-3,''event_chan'',1);';
% various spike counts from raster data
cmds{10} = 'spikecount_1400ms = raster2SpikeCount(raster,''start'',0,''end'',1.4);';
cmds{11} = 'spikecount_100ms = raster2SpikeCount(raster,''start'',0,''end'',0.1);';
cmds{12} = 'spikecount_50ms = raster2SpikeCount(raster,''start'',0,''end'',0.05);';
% smoothed receptive field
cmds{13} = 'smoothed_RF = smoother(spikecount_50ms);';
% stats vector for commands
stats = [(1:(length(cmds)))',zeros(length(cmds),1)];

dscrpt = '';


dets = struct('path',[pwd,'\'],'filename',{''},'modified',0,'backup',0,'date','','date_modified','');

%-------------------------------------------------------------
% create the SRC_ strings
%-------------------------------------------------------------
% Example:
%	if .src file name is "testsrcfile.src":
%		SRCrootname = 'SRC_testsrcfile'
%		SRCstructname = 'SRC_testsrcfile.'
%-------------------------------------------------------------
SRCrootname = ['SRC_' files.name(1:end-4)]
SRCstructname = ['SRC_' files.name(1:end-3)]

%-------------------------------------------------------------
% for current file, create empty data struct
%	data struct will have name "SRC_<filename without extension>
%-------------------------------------------------------------
evalstr(['SRC_', files.name(1:end-4),' = struct;'])

% create output filenames
% .src file
fname1 = [pwd filesep files.name(1:end-3) 'src'];
% sorted .plx file
fname2 = [pwd filesep files.name(1:end-4) '_plxsorted.plx'];
% efd file
fname3 = ['SRC_' files.name(1:end-3) 'efd'];

% date from file information
srcdate = [files.date];

% save src and plx file structs
[SRCstructname 'srcfiledetails = src;']
evalstr([SRCstructname 'srcfiledetails = src;']);
evalstr([SRCstructname 'plx_srcfiledetails = plx;'])

% save commands cell cmds in [SRCvarname.scriptcmds]
evalstr([SRCstructname 'scriptcmds = cmds;'])
evalstr([SRCstructname 'scriptstatus = stats;'])
evalstr([SRCstructname 'description = dscrpt;'])
evalstr([SRCstructname 'datafiledetails = dets;'])

evalstr([SRCstructname 'srcfiledetails.filename = {fname1};'])
evalstr([SRCstructname 'srcfiledetails.date = files.date;'])
evalstr([SRCstructname 'srcfiledetails.size = files.bytes;'])

evalstr([SRCstructname 'plx_srcfiledetails.filename = {fname2};'])

evalstr([SRCstructname 'datafiledetails.filename = fname3;'])
evalstr([SRCstructname 'datafiledetails.date = datestr(now);'])
evalstr([SRCstructname 'datafiledetails.date_modified = datestr(now);'])       


% save to .efd file
evalstr(['savePLXfile(' SRCrootname ',' '''' SRCstructname 'efd' '''' ');'])

% clear variable
evalstr(['clear ' SRCrootname])

