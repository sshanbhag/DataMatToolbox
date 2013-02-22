%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% initial setup
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% clear all vars and close all plots to avoid confusion
close all; clear all;
FORCEFLAG = 0;

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% file names
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

%------------------------------------------------------------
% set path and library depending on system
%------------------------------------------------------------
switch computer
	case {'PCWIN', 'PCWIN64'}
		DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
		datapath = 'F:\Work\Data\DataWave\DDFtest';
	case 'MACI64'
		datapath = '/Users/sshanbhag/Work/Data/DataWave/DDFtest';
end

%------------------------------------------------------------
% datafiles
%------------------------------------------------------------
% freq/response area, sorted spikes (bat)
filename = '829_01-05-2013--2729_FRA_Sorted.ddf';
% generate matfile name for output of converted .ddf data
[~, matfile] = fileparts(filename);
matfile = [matfile '.mat'];

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% convert DDF data to matfile (via NeuroShare) OR load matfile
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
if ~exist(fullfile(datapath, matfile), 'file') || FORCEFLAG
	fprintf('%s: converting ddf file to mat file...\n', mfilename)
	try
		D = DW.convertDDF2MAT(fullfile(datapath,filename), 'EVENT', 'SEGMENT', 'NEURAL');
	catch errMsg
		errMsg %#ok<NOPTS>
	end
	fprintf('%s: ....done\n', mfilename);
else
	fprintf('%s: loading .mat file %s\n', mfilename, fullfile(datapath, matfile));
	load(fullfile(datapath, matfile));
end

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% create Data Object from DDF neuroshare data matfile
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
d = DW.FRAdata(D, fullfile(datapath, filename));

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% plot FRA data
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
d.plotFRA(1, 255, [0 1000])

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%% save Data Object
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
save( fullfile(datapath, ['Dobj_' matfile]), 'd', 'D');
