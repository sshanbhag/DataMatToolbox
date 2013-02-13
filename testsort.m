%------------------------------------------------------------
%------------------------------------------------------------
%% initial setup
%------------------------------------------------------------
%------------------------------------------------------------
close all; clear all;
FORCEFLAG = 0;

%------------------------------------------------------------
%------------------------------------------------------------
%% file names
%------------------------------------------------------------
%------------------------------------------------------------

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
%********* BBN rate level function (bat)
% filename = '12-12-2012--2854_BBNrate.ddf';
%********* frequency tuning (bat)
% filename = '12-12-2012--2854_FreqScan2.ddf';
%********* vocal strings (bat)
% filename = '12-12-2012--2854_strings block.ddf';
%********* repetition rate (bat)
% filename = '12-12-2012--2854_RepRate0.ddf';
%********* syllables, sorted spikes (bat)
% filename = '01-03-2013--2961_syllable_block_new_sorted.ddf';
%********* syllables, sorted spikes, multispikes (bat)
% filename = '01-03-2013--2961_syllable_block_new_sorted.ddf';
%********* freq/response area, sorted spikes (bat)
% filename = '829_01-05-2013--2729_FRA_Sorted.ddf';
%********* repetition rate, sorted spikes (bat)
% filename = '829_01-05-2013--2729_repRate20_Sorted.ddf';
%********* Combination sensitivity, Multichannel recording (mouse)
% filename = '834_012913_combosens_1_spikes_sorted_1_small.ddf';
%********* multichannel (mouse)
filename = '834_020613_C_up_O_down_1_spikes_sorted.ddf'
% generate matfile name for output of converted .ddf data
[~, matfile] = fileparts(filename);
matfile = [matfile '.mat'];

%------------------------------------------------------------
%------------------------------------------------------------
%% convert to matfile or load matfile
%------------------------------------------------------------
%------------------------------------------------------------
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

%------------------------------------------------------------
%------------------------------------------------------------
%% create Data Object
%------------------------------------------------------------
%------------------------------------------------------------
d = DW.Data(D, fullfile(datapath, filename));

save( fullfile(datapath, ['Dobj_' matfile]), 'd', 'D');


