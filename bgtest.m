%------------------------------------------------------------
%------------------------------------------------------------
%% initial setup
%------------------------------------------------------------
%------------------------------------------------------------
% close figures 
close all; 
% clear variables and classes
clear classes;

% FORCEFLAG set to 1 will force conversion the .ddf file to a .mat file
% if the .mat file does not exist for the selected .ddf file, then the file will
% be automatically converted.
FORCEFLAG = 0;

% select the file to analyze from the datafiles list
FILENUM = 2;

% probenum and unitnum set to [] will cause script to ask for user input
% otherwise, set them to desired value
probenum = 1;
unitnum = 255;

% window for psth.  time in milliseconds relative to stimulus onset timestamp
psthwin = [-100 1000];

% background window to use (1 = beginning, 2 = ending)
bgwin_id = 1;

% duration of background window, in milliseconds
bgwin_duration = 100;

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
		% to use the DataWave Neuroshare library, set DLLName to the
		% path/filename of the DataWave Neuroshare .dll 
		% Note that the Neuroshare dll is only usable on Windows systems!
		DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
		% path to data files (Windows systems)
		datapath = 'F:\Work\Data\DataWave\DDFtest';
	case {'MACI64', 'GLNX86', 'GLNXA64'}
		% path to datafile for Mac or Linux
		datapath = '/Users/sshanbhag/Work/Data/DataWave/DDFtest';
	otherwise
		error('%s: unknown computer', mfilename);
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
	'827_12-18-2012--2649_syllable block_Sorted.ddf', ...
	'01-03-2013--2961_syllable_block_new_sorted.ddf', ...	% syllables, sorted spikes, multispikes (bat)
	'829_01-05-2013--2729_FRA_Sorted.ddf', ...		% freq/response area, sorted spikes (bat)
	'829_01-05-2013--2729_repRate20_Sorted.ddf', ...		% repetition rate, sorted spikes (bat)
	};
% filetype
%	RLF	= Rate Level Function (or any data that can fall under this
%										category)
%	FRA	= Frequency Response Area
battype = { ...
	'RLF', ...
	'RLF', ...
	'FRA', ...
	'RLF', ...
	};

% select the filename and datatype via the lists above and the FILENUM from
% inital config
filename = batfiles{FILENUM};
datatype = battype{FILENUM};

% generate matfile name for output of converted .ddf data
[~, matfile] = fileparts(filename);
matfile = [matfile '.mat'];

%------------------------------------------------------------
%------------------------------------------------------------
%% load Data Object
%------------------------------------------------------------
%------------------------------------------------------------
% write the data object.  for debugging use, the D struct 
% is also saved, but this will ultimately be redundant...
load( fullfile(datapath, ['Dobj_' matfile]));

%------------------------------------------------------------
%------------------------------------------------------------
%% grab background data
%------------------------------------------------------------
%------------------------------------------------------------
%----------------------------------------------------
% find the row of spiketimes that matches the unitnum
%----------------------------------------------------
unitrow = find(unitnum == d.Background(probenum).findClusters);
%----------------------------------------------------
% get the background spiketimes and waveforms for the desired 
% probe, unit and bg window (1 = first 30 seconds, 2 = final 30 seconds)
%----------------------------------------------------
% need to use cell2mat in order to get normal vector of timestamps 
% (instead of cell)
spikes = cell2mat(d.Background(probenum).t{unitrow, bgwin_id});
% get waveforms
waves = d.Background(probenum).wforms{unitrow, bgwin_id};

%----------------------------------------------------
% find random 100ms of bg spikes
%----------------------------------------------------
% generate a random start time from 0 to (30 seconds - bgwin_duration)
% note that the 30 seconds must be converted to milliseconds
% 
% the rand(n, m) function returns an [n X m] matrix of random numbers that are
% uniformly-distributed from 0 to 1.  To convert this value to one that covers
% the range of 0 - 30 seconds (minus window duration), simply scale by 
% the max value desired: in this case, 1000*30 - bgwin_duration.
tstart = (1000*30 - bgwin_duration) * rand(1, 1);
% convert to microseconds (timescale for timestamps)
tstart = 1000 * tstart;
% compute tend
tend = tstart + 1000*bgwin_duration;
% find spike timestamps within this window
[bgspikes, tmpindx] = find_valid_timestamps(spikes, tstart, tend, 0);
% get the valid waveforms
bgwaves = waves(cell2mat(tmpindx));



