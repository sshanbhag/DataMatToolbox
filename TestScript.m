%-----------------------------------------------------------------------------
% dwtester.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
%-----------------------------------------------------------------------------
% script to demonstrate ddf file/data access using Neuroshare toolbox and
% object-based analysis code in the DW package.
%
%-----------------------------------------------------------------------------
% See also: +DW package
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 25 February, 2013 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------


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
%% convert DDF data to matfile OR load matfile
%------------------------------------------------------------
%------------------------------------------------------------
if ~exist(fullfile(datapath, matfile), 'file') || FORCEFLAG
	% if the matfile for this ddf file does not exist, OR FORCEFLAG == 1, 
	% run convertDDF2MAT from the DW package to convert the .ddf data to a
	% Matlab struct (D) and save the data struct to a mat file via Neuroshare.  
	% The 'Event', 'Segment' and 'Neural' options tell the function to extract
	% the Event (markers), Segment (waveforms) and Neural (spike times)
	% Neuroshare entities
	fprintf('%s: converting ddf file to mat file...\n', mfilename)
	try
		D = DW.convertDDF2MAT(fullfile(datapath,filename), 'EVENT', 'SEGMENT', 'NEURAL');
	catch errMsg
		% trap any errors
		errMsg %#ok<NOPTS>
	end
	fprintf('%s: ....done\n', mfilename);
else
	% if the matfile exists, load it (saves time)
	fprintf('%s: loading .mat file %s\n', mfilename, fullfile(datapath, matfile));
	load(fullfile(datapath, matfile));
end

%------------------------------------------------------------
%------------------------------------------------------------
%% create Data Object
%------------------------------------------------------------
%------------------------------------------------------------
% check the datatype - if it is FRA (Frequency Response Area), use the 
% FRAdata class from the DW package to store/manage the data.  
% If the datatype is RLF, use the RateData class.  Initialize the classes 
% from the D struct obtained by convertDDF2MAT
if strcmpi(datatype, 'FRA')
	% load 
	d = DW.FRAdata(D, fullfile(datapath, filename));
	
elseif strcmpi(datatype, 'RLF')
	d = DW.RateData(D, fullfile(datapath, filename));

else
	error('%s: Unknown data type %s', mfilename, datatype);
end

%------------------------------------------------------------
%------------------------------------------------------------
%% plot waveforms
%------------------------------------------------------------
%------------------------------------------------------------
if isempty(probenum) || isempty(unitnum)
	% plot unit waveforms (overlaid) for all probes and units
	d.plotUnitWaveforms;
else
	% plot unit waveforms (overlaid) for selected probe and unit
	d.plotUnitWaveforms('probe', probenum, 'unit', unitnum);
end


%------------------------------------------------------------
%------------------------------------------------------------
%% get probe and unit number from user
%------------------------------------------------------------
%------------------------------------------------------------
% Select Probe if probenum is empty
if isempty(probenum)
	qFlag = 1;
	while qFlag
		fprintf('Probes found:\t')
		fprintf('%d ', 1:d.Nprobes);
		fprintf('\n');
		probenum = input('Select Probe: ');
		if isempty(probenum)
			fprintf('must select a probe number!\n\n')
		elseif ~any(probenum == 1:d.Nprobes)
			fprintf('please select a probe in range!\n\n')
		else
			qFlag = 0;
		end
	end
end	% END if 
% Select Unit
if isempty(unitnum)
	qFlag = 1;
	while qFlag
		fprintf('\n\nUnits (clusters) for Probe %d:\t', probenum);
		fprintf('%d ', d.Probes(probenum).cluster);
		fprintf('\n');
		unitnum = input('Select Unit (cluster #): ');
		if isempty(unitnum)
			fprintf('must select a unit number!\n\n')
		elseif ~any(unitnum == d.Probes(probenum).cluster)
			fprintf('please select a unit that is in range!\n\n')
		else
			qFlag = 0;
		end
	end
end	% END if

%------------------------------------------------------------
%------------------------------------------------------------
%% plot appropriate to type
%------------------------------------------------------------
%------------------------------------------------------------
if strcmpi(datatype, 'FRA')
	% plot FRA data, using window from 0 to 1000 msec
	% 0 corresponds to the stimulus onset marker/timestamp
	d.plotFRA(probenum, unitnum, [0 1000]);

elseif strcmp(datatype, 'RLF')
	% plot rasters/psth for each stimulus and attenuation
	% individual stimuli will be plotted in separate figures (for all atten
	% values)
	d.plotRasterAndPSTH('probe', probenum, 'unit', unitnum, 'offset', [-100 0])
% 	spikecount = d.countSpikes(probenum, unitnum, [0 50]);
else
	error('%s: Unknown data type %s', mfilename, datatype);
end

%------------------------------------------------------------
%------------------------------------------------------------
%% get 1 ms psths 
%------------------------------------------------------------
% Format:
% 	PSTH:	{# atten vals, # stim groups}, where each element
% 			is a matrix of size [# sweeps (reps), # bins]
% 
% 	bins:	[1, length(psthwin(1):binsize:psthwin(2))] vector of 
% 				time bins for PSTH.
% 
% 	spiketimes:	{# attenvals, # stim groups}  cell matrix, each element is
% 					a cell vector, {# sweeps}, each containing a vector of
% 					spiketimes (in milliseconds) for that sweep;
% 
%	StimInfo		struct
% 		StimInfo.AttenVals	{# attenvals, # stim groups} array
% 									of attenuation values for each atten/stim group
% 									combination
% 		StimInfo.VarValues	{# attenvals, # stim groups} array
% 									of stimulus variable.  will depend on stimulus
% 									type.  Wav stimuli will have .wav filename, 
% 									Tone will have Frequency, noise will have
% 									low and upper freq vals
%------------------------------------------------------------
%------------------------------------------------------------
binsize = 1;	% msec
[PSTH, bins, spiketimes, StimInfo] = d.computePSTH(probenum, unitnum, binsize, psthwin);


%------------------------------------------------------------
%------------------------------------------------------------
%% save Data Object
%------------------------------------------------------------
%------------------------------------------------------------
% write the data object.  for debugging use, the D struct 
% is also saved, but this will ultimately be redundant...
save( fullfile(datapath, ['Dobj_' matfile]), 'd', 'D');


