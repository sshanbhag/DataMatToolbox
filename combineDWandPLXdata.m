%-----------------------------------------------------------------------------
% combineDWandPLXdata
%-----------------------------------------------------------------------------
% ExpFun Project
%-----------------------------------------------------------------------------
%
% Combines DataWave text files (stimulus info) and Plexon timestamp data
% from .plx files for subsequent use in ExpFunBrowser program
% 
%-----------------------------------------------------------------------------
% Structures created by combineDWandPLXdata:
% 
% 	stimdata		Nx1 struct array of stimulus information
% 		Field						Description
%		 name						 .wav file name for stimulus, or frequency for freq-resp area
% 		 indices					 locations in VARS arrays for this stimulus
% 		 nsweeps					 # of sweeps for this stimulus
% 		 tstamps					 timestamps for this stimulus
% 		 attenvals				 attenuation values used with this stimulus
% 		 n_attenvals			 # of attenuation values used with this stimulus
% 		 sweepindex				 (n_attenval X 1) cell array containing indices for 
% 		 							  all attenuation values for this stimulus. 
% 		 sweeps_per_atten		 # sweeps for each attenuation level (for this stimulus)
% 		 tstamps_per_atten	 # of timestamps  for each atten and stimulus combination
% 		 sweepstart_t			 times (in seconds) for the start of each sweep
% 		 sweepend_t				 times (seconds) for the end of each sweep
% 		 spikes					 spike timestamps (seconds) for each stimulus and 
% 		 							  attenuation combination
% 		 stim_idnum				 Stimuli are assigned an id number because some of
% 									  expfunbrowser's functions barf on character-valued
% 									  stimulus id's.  This stimulus' id value is stored here
% 
% 
%	vars		structure with information about stimulus variables
%  	Field				Description
% 		 values			 [N X m] vector of marker tags
% 		 names			 {1 X m} cell vector of text description of values
% 							  (e.g., {'atten'  'wavfilename'}
% 		 description	 text description of data
% 
% 	plx		Plexon timestamp data structure 
% 	         **note that event information is more than likely to be useless**
%  	Field						Description
% 		 filename				 string containing full path and name of .plx file
% 		 ts						 {1 X total_noof_units} cell array of timestamps in seconds
% 		 total_noof_units		 # of units sorted
% 		 unit_id					 
% 		 chan_id
% 		 ts_counts				 total number of spikes per unit
% 		 event_ts				 event timestamps
% 		 noof_event_chans
% 		 event_chan_id
% 		 event_ts_counts 
% 
% 	rawdata		structure containing raw information from the DataWave text file
% 		Field			Description
% 		 dwdata		 datawave event data
% 		 marker		 datawave marker data in vector form
% 
% 	spikes			expfun browser spike structure
% 	 Example:
% 		values: [16x5x3 double]
% 		interval_count: [16x5x3 double]
% 		sd: [16x5x3 double]
% 		conf90: [16x5x3 double]
% 		conf95: [16x5x3 double]
% 		dimvals: {[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]  [10 30 50 70 80]  [1 2 3]}
% 		noof_dims: 3
% 		dim_names: {'Unit No.'  'Attenuation'  'Filename'}
% 		total_noof_units: 16
% 		unit_id: [16x1 double]
% 		chan_id: [16x1 double]
% 		plot_mode: 'spike_count'
% 		vars_description: 'Test Stimuli'
% 		plx_files: '/Data/DataWave/TestData/EXPFUNTESTER/EXPFUN_test_1per2sec_cahn1to7_sorted.plx'
% 		fn_args: {'var_matrix'  'variables'  [1x1 struct]  'end'  [1]  'start'  [0]}
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% See Also: expfunbrowser
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Code History
%-----------------------------------------------------------------------------
% Revisions:
%	17 May, 2011 (SJS): file created
%-----------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% first, get file information from user
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
[dwfilename, dwpath] = uigetfile('*.txt', 'Select DataWave txt file...')

if isequal(dwfilename, 0) || isequal(dwpath, 0)
	% user selected 'CANCEL' so end
	disp('Cancelling...')
	return
end

% set the putative plxpath to the dwpath (assume that the plexon file is in the
% same directory as the datawave text file)
plxpath = dwpath;
% build the plx filename
[tmp, filename, fext] = fileparts(dwfilename);
plxfilename = [filename '.plx'];

if ~exist(fullfile(plxpath, plxfilename), 'file')
	% plexon file not found
	warning('%s: could not locate matching plx file %s in directory %s !!!!!!!', ...
					mfilename, plxfilename, plxpath);
	% give user a chance to manually locate the plexon file
	[plxfilename, plxpath] = uigetfile('*.plx', ...
										['Select Plexon file for DataWave file ' dwfilename ' ...']);

	if isequal(plxfilename, 0) || isequal(plxpath, 0)
		% user selected 'CANCEL' so end
		disp('Cancelling...')
		return
	end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% initialize plx file details struct (presently, this is unused but left
% in this code for future use
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
plxfiledetails = struct(	'noof_files',			1, ...
										'filename',				{{plxfilename}}, ...
										'startinterval',		1, ...
										'noof_intervals',		-1	);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% load the data using loadDWStimData function
%--------------------------------------------------------------------------
[stimdata, vars, plx, rawdata] = loadDWStimData(fullfile(plxpath, plxfilename), fullfile(dwpath, dwfilename));

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% convert data spike timestamps into spikes struct format used by
% expfunbrowser.
% Technique suggested by Chris Sumner
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% need to set sweep duration
%--------------------------------------------------------------------------
loopFlag = 1;
while loopFlag
	tmp = input('Enter sweep duration in seconds: ', 's');
	if isempty(tmp)
		loopFlag = 1;
	else
		tmpnum = str2num(tmp);
		if isempty(tmpnum)
			disp(['invalid value for sweep duration (' tmp ')!'])
			loopFlag = 1;
		elseif tmpnum <= 0
			disp(['Error: sweep duration must be greater than 0'])
			loopFlag = 1;
		else
			loopFlag = 0;
		end
	end
end
sweepDuration = tmpnum;

spikes = PLX2SpikeCount(plx,'var_matrix','variables',vars,'end', sweepDuration);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plot data using expfunbrowser command displayGenData
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
display_GenData(spikes, 'xdim',2, 'axesdimlist',1, 'minval',10,'title', vars.description);

