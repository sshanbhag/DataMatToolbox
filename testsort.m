%------------------------------------------------------------
%------------------------------------------------------------
%% initial setup
%------------------------------------------------------------
%------------------------------------------------------------
close all; clear all;
FORCEFLAG = 0;
FILENUM = 6;
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

% MOUSE
%********* Combination sensitivity, Multichannel recording (mouse)
%**** CRASHES NEUROSHARE MATLAB API (not NeuroExplorer)
% filename = '834_012913_combosens_1_spikes_sorted_1_small.ddf'
%********* multichannel (mouse)
% filename = '834_020613_C_up_O_down_1_spikes_sorted.ddf';

% BAT
batfiles = { ...
	'12-12-2012--2854_BBNrate.ddf', ...		%	 BBN rate level function
	'12-12-2012--2854_FreqScan2.ddf', ...	% frequency tuning (bat)
	'12-12-2012--2854_strings block.ddf', ...		% vocal strings (bat)
	'12-12-2012--2854_RepRate0.ddf', ...		% repetition rate (bat)
	'01-03-2013--2961_syllable_block_new_sorted.ddf', ...	% syllables, sorted spikes, multispikes (bat)
	'829_01-05-2013--2729_FRA_Sorted.ddf', ...		% freq/response area, sorted spikes (bat)
	'829_01-05-2013--2729_repRate20_Sorted.ddf', ...		% repetition rate, sorted spikes (bat)
	};
filename = batfiles{FILENUM};

% generate matfile name for output of converted .ddf data
[~, matfile] = fileparts(filename);
matfile = [matfile '.mat'];

%------------------------------------------------------------
%------------------------------------------------------------
%% convert DDF data to matfile OR load matfile
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
d = DW.FRAdata(D, fullfile(datapath, filename));


%------------------------------------------------------------
%------------------------------------------------------------
%% save Data Object
%------------------------------------------------------------
%------------------------------------------------------------
save( fullfile(datapath, ['Dobj_' matfile]), 'd', 'D');


