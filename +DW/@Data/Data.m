%-----------------------------------------------------------------------------
% Data.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW Package
% Class Definition
%-----------------------------------------------------------------------------
%	Data facilitates access to data output by loadDWfile.m
%
%  obj = Data('xxx_y_zz_ttt_n.mat') loads the following:
%     obj.fpath;   %path to file
%     obj.fname;   %name of file
%     obj.fext;    %file extension
%     obj.Background; %really this should be a class too...
%     obj.D;          %really this should be a class too...
%     obj.Stimulus;   %really this should be a class too...
%     obj.fullfname;  %full file name
%-----------------------------------------------------------------------------
% See also: loadDWfile (function)
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% initial coding and design:
%	Tony Slagle
%	tonyslagle@gmail.com
%
% Continuing development: 
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: January, 2012 (TS)
%
% Revisions:
%	17 May, 2012 (SJS): 
%	 -	updated documentation
%	 -	renamed file and class from LoadDWfileData.m to Data.m
%	 -	beginning massive Objectification...
%	13 Dec 2012 (SJS):
%	 -	fixing issues on initialization of object
%	14 Jan 2013 (SJS)
%	 -	major overhaul to use NeuroShare Matlab API and DataWave Neuroshare
%		interface object (@NS)
%	31 Jan 2013 (SJS): reworking to avoid segmentation faults 
%	12 Feb 2013 (SJS): added plotting routines
%	22 Feb 2013 (SJS): modifying get spikes to allow for pre/post spike time
%-----------------------------------------------------------------------------
% TO DO: lots....
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% classdef (ConstructOnLoad = true) Data < handle
classdef Data < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		fpath			%path to file
		fname			%name of file
		fext			%file extension
		Info			% DWinfo object
		Markers		% Marker object array
		Background	% really this should be a class too...
		Stimuli		% Stimulus object array
		Probes		% Probe object array
		Units			% Unit object array
		Nmarkers
		Nprobes
		Nunits
	end
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define other properties
	%------------------------------------------------------------------------
	properties
		fullfname;  %full file name
	end
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
  
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		% Initializes the object
		%------------------------------------------------------------------------
		function obj = Data(varargin)
		%---------------------------------------------------------------------	
		%	Data(<fileName>) 
		%	Constructor method
		%	opens file called fileName (char) or opens
		%  a dialog box to get a filename if the fileName provided does not exist.
		%---------------------------------------------------------------------	

			%--------------------------------------------------------
			% Check inputs
			%--------------------------------------------------------
			% if no inputs, return
			if isempty(varargin)
				return
			end
			
			% if 1 input, need to perform a few checks
			if length(varargin) == 1
				% check if varargin{1} was directly passed in from
				% a subclass (using obj@DW.Stimulus(varargin) syntax)
				if strcmpi(inputname(1), 'varargin')
					% make local copy of varargin{1}
					tmp = varargin{1};
					if length(tmp) == 1
						% tmp is D structure
						obj.initFromStruct(tmp{1});
					elseif length(tmp) == 2
						% caller provided D structure and filename
						obj.initFromStruct(tmp{1});
						obj.fullfname = tmp{2};
						[obj.fpath, obj.fname, obj.fext] = fileparts(tmp{2});
					else
						fprintf('%s: strange inputs... ', mfilename);
						fprintf('%s\t', tmp);
						fprintf('\n');
						error('%s: invalid inputs', mfilename);
					end
					clear tmp
				else
					% only D struct provided
					obj.initFromStruct(varargin{1});
				end
			elseif length(varargin) == 2
				% caller provided D struct and filename
				obj.initFromStruct(varargin{1});
				obj.fullfname = varargin{2};
				[obj.fpath, obj.fname, obj.fext] = fileparts(varargin{2});
			else
				error('%s: invalid input args', mfilename);
			end
		end	% END Data CONSTRUCTOR
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Data input/reading methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function initFromStruct(obj, D)
		%------------------------------------------------------------------------
		% Data.initFromStruct
		%------------------------------------------------------------------------
		%	Takes information from NeuroShare data struct D and initializes
		%	Data object properties 
		%------------------------------------------------------------------------
		
			%-------------------------------------------------------
			% store info bits in obj.Info!
			%-------------------------------------------------------
			infotags = {'FileType', 'EntityCount', 'TimeStampResolution', ...
								'TimeSpan', 'AppName', 'Time_Year', 'Time_Month', ...
								'Time_Day', 'Time_Hour', 'Time_Min', 'Time_Sec', ...
								'Time_MilliSec', 'FileComment' };
			for n = 1:length(infotags)
				obj.Info.(infotags{n}) = D.(infotags{n});
			end
			%-------------------------------------------------------
			% store Markers
			%-------------------------------------------------------
			if ~isfield(D, 'Event')
				fprintf('%s Error: no Event entities in %s\n', ...
																		mfilename, inputname(1));
				return;
			elseif D.nEvent == 0
				fprintf('%s Error: nEvent = %d\n', mfilename, D.nEvent);
				return;
			else
				obj.loadMarkers(D.Event);
			end
			%-------------------------------------------------------
			% convert to stimuli
			%-------------------------------------------------------
			obj.loadStimuli;
			%-------------------------------------------------------
			% get sweep times
			%-------------------------------------------------------
			obj.Stimuli.extractTimeFromMarkers(obj.Markers);
			%-------------------------------------------------------
			% find stimulus groups (same apart from attenuation)
			%-------------------------------------------------------
			obj.Stimuli.findCommon;
			%-------------------------------------------------------
			% convert segments to Probes
			%-------------------------------------------------------
			if ~isfield(D, 'Segment')
				fprintf('%s Error: no Segment entities in %s\n', ...
																		mfilename, inputname(1));
				return;
			elseif D.nSegment == 0
				fprintf('%s Error: nSegment = %d\n', mfilename, D.nSegment);
				return;
			else
				obj.loadProbesFromSegment(D.Segment);
				obj.loadBackground;
			end
		end	% END initFromStruct
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function varargout = loadMarkers(obj, Events)
		%------------------------------------------------------------------------
		% [Events, errFlg] = Data.loadMarkers
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% Initial setup and checks
			%-----------------------------------------------------------
			DW.DataWaveDefaults	% load defaults
			errFlg = 0;
			if isempty(Events)
				errFlg = 1;
				warning('%s: no Events!', mfilename);
				return
			else
				nEvents = length(Events);
			end

			%-----------------------------------------------------------
			% parse Markers
			%-----------------------------------------------------------			
			% From Neuroshare/DataWave, each channel's stimulus Marker
			% is stored in an individual Event entity. the left (channel
			% 1) marker is in EventID = 1, the right (channel 2) is in 
			% EventId 2.
			% 
			% create a list of Event IDS for each event entity.
			% 
			%	EventID == 0 --> non-marker entity
			%	EventID == 1 --> L channel markers
			%	EventID == 2 --> R channel markers
			%
			% 7 Mar 2013 (SJS): issue: some inconsistency in files.
			%	need to rethink assignment, since some files do not
			%	have L (or R) markers, and thus only 2 Events.
			%	
			%	Also, the # of L and R markers may not match!
			%
			%-----------------------------------------------------------			
			EventID = zeros(nEvents, 1);

			% first, find the R and L marker events (R is usually Events(2), 
			% L is usually Events(3), but try not to assume)
			for v = 1:nEvents
				% scan the CSVDesc from the NeuroShare Events.Info to get
				% a list (usually partial - for sound stimulus markers, the 
				% DataWave Neuroshare (or Matlab API, unknown at moment...) leaves
				% out some of the labels
				tmp = textscan(Events(v).Info.CSVDesc, '%s', ...
												'Delimiter', ',', ...
												'MultipleDelimsAsOne', 1);
				field_names = tmp{1};
				clear tmp;
				% check if the current Event entity has 'SoundType' as part of
				% the first element of its fields
				if strncmpi(field_names{1}, 'SoundType', length('SoundType'))
					% if so, assume that this is a "Marker" Entity and check for
					% R or L to determine Right or Left sound marker
					if strcmpi(field_names{1}(end), 'L')
						% left marker
						EventID(v) = L;
					elseif strcmpi(field_names{1}(end), 'R')
						% right marker
						EventID(v) = R;
					else
						EventID(v) = 0;
					end
				end
			end

			evL = find(EventID == L);
			evR = find(EventID == R);
			
			fprintf('Left Markers -> Event(%d)\n', evL);
			fprintf('Right Markers -> Event(%d)\n', evR);
			
			% check to ensure that both L and R Entities were found
			if ~all( any(EventID == L) & any(EventID == R) )
				% if not abort
				error('%s: Markers for both (L and R) channels not found!', ...
																						mfilename);
			end
			
			% check MarkerEvents
			MarkerEvents = find(EventID);
			if length(MarkerEvents) ~= 2
				error('%s: # of MarkerEvents is %d (must be equal to 2)',...
								mfilename, length(MarkerEvents));
			end
			evCount = [Events(MarkerEvents(1)).EventCount ...
										Events(MarkerEvents(2)).EventCount];
			% make sure # of events are the same
			if ~(evCount(1) == evCount(2))
				warning('%s: EventCount mismatch between L and R channels!', ...
																						mfilename);
				[EventCount, minindx] = min(evCount);
				obj.Nmarkers = EventCount;
				fprintf('Using lower count (%d)\n', obj.Nmarkers);
			else
				EventCount = Events(MarkerEvents(1)).EventCount;
				obj.Nmarkers = EventCount;
			end



			%-----------------------------------------------------------			
			% allocate Markers object array - need to figure out soln
			% for only 1 channel of markers....
			%-----------------------------------------------------------			
			% Matlab thinks obj.Markers is a double, and will throw an
			% an error.  so, we trick Matlab in order to do so:  first, need to
			% initialize Markers as a single marker, then allocate 
			% full array. 
			obj.Markers = DW.Marker;
			obj.Markers(EventCount, 1) = DW.Marker;
			% loop through the events
			for n = 1:EventCount
				% parse strings into cell arrays; to preserve empty fields, do NOT
				% treat successive delimiters/whitespace as one (in call to csvscan)
				tmpR = csvscan(Events(evR).Data{n}, 0);
				tmpL = csvscan(Events(evL).Data{n}, 0);
				% ASSUME that one of the lists will not have an outputfile field.
				% pad with empty values 
				if length(tmpR) ~= MARKER_NBASE
					dlen = MARKER_NBASE - length(tmpR);
					if dlen < 0
						error('%d: length(tmpR) > MARKER_NBASE');
					else
						tmpadd = cell(dlen, 1);
						for t = 1:dlen
							tmpadd{t} = '';
						end
						tmpR = [tmpR; tmpadd]; %#ok<AGROW>
						clear tmpadd
					end
				end
				if length(tmpL) ~= MARKER_NBASE
					dlen = MARKER_NBASE - length(tmpL);
					if dlen < 0
						error('%d: length(tmpL) > MARKER_NBASE');
					else
						tmpadd = cell(dlen, 1);
						for t = 1:dlen
							tmpadd{t} = '';
						end
						tmpL = [tmpL; tmpadd]; %#ok<AGROW>
						clear tmpadd
					end
				end
				% create a single cell array of strings
				tmp = [tmpR; tmpL];
				% convert tags appropriately, first allocating a list for values
				elist = cell(size(tmp));
				% loop through the marker tags
				for t = 1:MARKER_NMARKERS
					% check for type
					switch MARKER_TYPES{t} %#ok<USENS>
						case 'char'
							elist{t} = tmp{t};
						case 'float'
							elist{t} = str2double(tmp{t});
						case 'int'
							elist{t} = str2double(tmp{t});
						otherwise
							elist{t} = tmp{t};
					end
				end	% END t -> MARKERS_NMARKERS

				% assign values to marker object
				obj.Markers(n).setValuesFromEventList(elist);
				% set timestamp
				obj.Markers(n).setTimestamp([	Events(evL).TimeStamp(n)	...
														Events(evR).TimeStamp(n) ]);
				% set ID
				obj.Markers(n).setID(n);
				clear elist
			end
			clear tmpR tmpL tmp elist
			%-------------------------------
			% assign outputs
			%-------------------------------
			if nargout > 0
				if nargout <= 1
					varargout{1} = errFlg;
				end
				if nargout > 1
					varargout{2} = Events;
				end
			end
		end	%END loadMarkers
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		function varargout = loadStimuli(obj)
		%------------------------------------------------------------------------
		% Stimuli = Data.loadStimuli
		%------------------------------------------------------------------------
			DW.DataWaveDefaults	% load defaults
			%-----------------------------------------------------------
			% load markers if they haven't been loaded yet
			%-----------------------------------------------------------
			if isempty(obj.Markers)
				obj.loadMarkers
			end
			if ~obj.Nmarkers
				error('%s: no markers in Data!')
			end
			%-----------------------------------------------------------
			% Initialize and build the StimulusList object (obj.Stimuli)
			% from the Markers
			%-----------------------------------------------------------
			obj.Stimuli = DW.StimulusList(obj.Markers);
			%-----------------------------------------------------------
			% assign outputs
			%-----------------------------------------------------------
			if nargout == 0
				return
			else
				varargout{1} = obj.Stimuli;
			end
		end	% END loadStimuli
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function varargout = loadProbesFromSegment(obj, Segment)
		%------------------------------------------------------------------------
		% Probes = Data.loadProbesFromSegment
		%------------------------------------------------------------------------
			DW.DataWaveDefaults	% load defaults
			%-----------------------------------------------------------
			% load markers if they haven't been loaded yet
			%-----------------------------------------------------------
			if isempty(obj.Markers)
				obj.loadMarkers
			end
			if ~obj.Nmarkers
				error('%s: no markers in Data!')
			end
			% check if there are any segment entities in the data
			if isempty(Segment)
				warning('%s: there are no Segment entities in data!', mfilename);
				return
			else
				fprintf('%s: Loading probe data from Segment entities\n', ...
									mfilename);
				nSegments = length(Segment);
			end
		
			%-----------------------------------------------------------
			% the segment entities have the spiketimes
			% as well as the spike snippets (if stored).  The Neural
			% entitites could be used, but are somewhat redundant and 
			% do not have the snippet waveforms.  
			%-----------------------------------------------------------
			obj.Nprobes = nSegments;			
			% create the Probe object
			obj.Probes = DW.Probe;
			if nSegments > 1
				obj.Probes(nSegments, 1) = DW.Probe;
			end	
			
			% pull out timestamps based on unit ID
			% according to Neuroshare, unit 0 is uncategorized, 255 is noise
			% this is not necessarily true
			for n = 1:obj.Nprobes
				% store timestamps for units with ID 0 - should
				% probably figure out a way to store/sort all unitIDs at some
				% point........
				% also, convert to microseconds (from seconds!!!)
				obj.Probes(n).name = Segment(n).SourceInfo.ProbeInfo;
				obj.Probes(n).samprate = Segment(n).Info.SampleRate;
				obj.Probes(n).cluster = unique(Segment(n).UnitID);
				obj.Probes(n).nclusters = length(obj.Probes(n).cluster);
				obj.Probes(n).t = cell(obj.Probes(n).nclusters, 1);
				obj.Probes(n).wforms = cell(obj.Probes(n).nclusters, 1);
				obj.Probes(n).time_units = 1e-6;
				obj.Probes(n).nwindows = 1;
				% loop through units - id 0 will usually be first, junk is 255
				for u = 1:obj.Probes(n).nclusters
					% convert timestamp to usec from seconds
					obj.Probes(n).t{u} = ...
									1e6 * Segment(n).TimeStamp(...
										Segment(n).UnitID == obj.Probes(n).cluster(u));
					obj.Probes(n).wforms{u} = ...
									Segment(n).WaveForm(...
										Segment(n).UnitID == obj.Probes(n).cluster(u));
				end	% END Nclusters loop
			end	% END nSegments loop
			
			clear Segment
			if nargout
				varargout{1} = unitID;
			end
		end	% END loadProbesFromSegmentData
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		function varargout = loadProbesFromNeural(obj, Neural)
		%------------------------------------------------------------------------
		% Probes = Data.loadProbesFromNeural
		%	uses neural entities to load probe information.  this is less
		% preferable than getting unit info from the Segment entities, since
		% the Neural entities do not have waveform snippets that can be used
		% to examine the unit quality!
		%------------------------------------------------------------------------
			DW.DataWaveDefaults	% load defaults
			%-----------------------------------------------------------
			% load markers if they haven't been loaded yet
			%-----------------------------------------------------------
			if isempty(obj.Markers)
				obj.loadMarkers
			end
			if ~obj.Nmarkers
				error('%s: no markers in Data!')
			end
			% check if there are any segment entities in the data
			if isempty(Neural)
				warning('%s: there are no Neural entities in data!', mfilename);
				return
			else
				nNeural = length(Neural);
			end
		
			%-----------------------------------------------------------
			% you'd assume that the neural events from 
			% DataWave/Neuroshare are what we want for the spike data.
			% you'd be wrong.  the segment entities have the spiketimes
			% as well as the spike snippets (if stored).
			%-----------------------------------------------------------
			% create the Probe object
			obj.Probes = DW.Probe;
			if nNeural > 1
				obj.Probes(nNeural, 1) = DW.Probe;
			end	
			
			% pull out timestamps based on unit ID
			% according to Neuroshare, unit 0 is uncatergorized, 255 is noise
			for n = 1:nNeural
				% store timestamps for units with ID 0 - should
				% probably figure out a way to store/sort all unitIDs at some
				% point........
				% also, convert to microseconds (from seconds!!!)
				keyboard
				obj.Probes(n).name = Neural(n).Info.ProbeInfo;
				obj.Probes(n).cluster = Neural(n).Info.SourceUnitID;
				obj.Probes(n).nclusters = 1;
				obj.Probes(n).t = Neural(n).TimeStamps;
				obj.Probes(n).time_units = 1e-6;
				obj.Probes(n).nwindows = 1;
			end	% END nSegments loop
			obj.Nprobes = nNeural;
			if nargout
				varargout{1} = obj.Probes;
			end
		end	% END loadProbesFromNeuralData
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function varargout = loadBackground(obj, varargin)
		%------------------------------------------------------------------------
		% Probes = Data.loadBackground
		%------------------------------------------------------------------------
			fprintf('%s: loading Background spikes from Probes...\n', mfilename);
			%-----------------------------------------------------------
			% make sure probes are loaded!
			%-----------------------------------------------------------
			if isempty(obj.Nprobes)
				error('%s: no probes loaded!', mfilename);
			end
			% default bg windows
			bgwindow =	{	[0 30], ...
								[(obj.Info.TimeSpan - 30)	obj.Info.TimeSpan], ...
							};
			% if user provided them, use them
			if nargin == 2
				bgwindow = varargin{1};
			end
			nwindows = length(bgwindow);
			
			%-----------------------------------------------------------
			% create Background object
			%-----------------------------------------------------------
			obj.Background = DW.Probe;
			if obj.Nprobes > 1
				obj.Background(obj.Nprobes, 1) = DW.Probe;
			end
			
			% more or less, copy probe data.  need to get timestamps for each
			% unit for the background windows
			for n = 1:obj.Nprobes
				% also, convert to microseconds (from seconds!!!)
				obj.Background(n).name = obj.Probes(n).name;
				obj.Background(n).samprate = obj.Probes(n).samprate;
				obj.Background(n).cluster = obj.Probes(n).cluster;
				obj.Background(n).nclusters = obj.Probes(n).nclusters;
				obj.Background(n).time_units = 1e-6;
				obj.Background(n).nwindows = nwindows;
				
				% find spikes within windows
				obj.Background(n).t = cell(obj.Probes(n).nclusters, nwindows);
				obj.Background(n).wforms = cell(obj.Probes(n).nclusters, nwindows);
				
				% loop through units - id 0 will usually be first, junk is 255
				for u = 1:obj.Probes(n).nclusters
					for w = 1:nwindows
					% convert timestamp to usec from seconds
						[obj.Background(n).t{u, w}, tempi] = ...
								find_valid_timestamps(	obj.Probes(n).t{u}, ...
																1e6 * bgwindow{w}(1), ...
																1e6 * bgwindow{w}(2)	);
						obj.Background(n).wforms{u, w} = ...
														obj.Probes(n).wforms{u}(tempi{1});
					end	% END nwindows loop
				end	% END Nclusters loop
			end	% END nSegments loop
			if nargout
				varargout{1} = obj.Background;
			end
					
		end	% END function
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Analysis methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		function varargout = computePSTH(obj, probenum, unitnum, binsize, psthwin)
		%------------------------------------------------------------------------
		% computes psth for spike data
		%------------------------------------------------------------------------
		% H = Data.computePSTH(obj, probenum, unitnum, binsize, psthwin)
		%	
		%	If spikes are desired from before the stim onset timestamp or after 
		%	the next stim onset timestamp, an 'offset' option of form
		%	[pretime posttime] may be included.  times must be in milliseconds
		%------------------------------------------------------------------------

			%------------------------------------------------
			% ensure that inputs are in bounds
			%------------------------------------------------
			if ~between(probenum, 1, obj.Nprobes)
				error('%s: probe must be in range [1:%d]', mfilename, obj.Nprobes);
			end
			if ~any(unitnum == obj.Probes(probenum).cluster)
				error('%s: unit not found in Probe]', mfilename, probenum);
			end
			if binsize <= 0
				error('%s: binsize must be greater than 0 ms', mfilename)
			end
			fprintf('%s: computing psth for probe %d, unit %d\n', ...
												mfilename, probenum, unitnum);
			%--------------------------------------------------
			% # of groups available
			%--------------------------------------------------
			ngroups = length(obj.Stimuli.GroupList);
			%--------------------------------------------------
			% # of atten available
			%--------------------------------------------------
			natten = zeros(ngroups, 1);
			for g = 1:ngroups
				Sindx = obj.Stimuli.GroupList{g};
				natten(g) = length(Sindx);
			end
			maxnatten = max(natten);
			%--------------------------------------------------
			% allocate arrays
			%--------------------------------------------------
			% compute # bins based on binsize and psthwin
			bins = psthwin(1):binsize:psthwin(2);
			nbins = length(bins);
			% psth data in H will be {maxnatten, ngroups} with each element
			% being a [nreps X nbins] array
			H = cell(maxnatten, ngroups);
			allspikes = cell(maxnatten, ngroups);
			%------------------------------------------------
			% loop through groups
			%------------------------------------------------
			% loop through # of stimulus groups (filenames, frequencies, etc.)
			for g = 1:ngroups
				tmpspikes = obj.getSpikesForStim(g, ...
															probenum, unitnum, ...
															'window', psthwin);
				% loop through # of attenuation values for this stim
				for n = 1:natten(g)
					% get spikes for this stim/atten level
					allspikes{n, g} = tmpspikes{n};
					nsweeps = length(allspikes{n, g});
					H{n, g} = zeros(nsweeps, nbins);
					% loop through sweeps
					for t = 1:nsweeps
						% convert spiketimes to milliseconds 
						%  (required for psth() func)
						allspikes{n, g}{t} = 0.001*allspikes{n, g}{t};
						% compute psth
						[H{n, g}(t, :), tbins] = psth(allspikes{n, g}{t}, binsize, psthwin);						
					end
				end
			end	% END g
			if any(nargout == 1:4)
				varargout{1} = H;
			end
			if any(nargout == 2:4)
				varargout{2} = bins;
			end
			if any(nargout == 3:4)
				varargout{3} = allspikes;
			end
			if nargout == 4
				varargout{4} = obj.getStimInfo;
			end
		end	% END computePSTH		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
	
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Information methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
			
		%------------------------------------------------------------------------
		function StimInfo = getStimInfo(obj)
		%------------------------------------------------------------------------
		
			% get # groups
			ngroups = length(obj.Stimuli.GroupList);
			% assume same # of levels for all groups!
			nlevels = length(obj.Stimuli.GroupList{1});
			% allocate StimInfo bits
			StimInfo.VarValues = cell(nlevels, ngroups);
			StimInfo.AttenValues = cell(nlevels, ngroups);
			
			% loop through groups
			for g = 1:ngroups
				nlevels = length(obj.Stimuli.GroupList{g});
				% loop through levels
				for n = 1:nlevels
					% get list of Stimuli for this group
					gIndx = obj.Stimuli.GroupList{g}(n);
					% get channel index (1 = L, 2 = R)
					cIndx = DW.chan2indx(obj.Stimuli.Channel(gIndx));
					% get Attenuation value
					StimInfo.AttenValues{n, g} = obj.Stimuli.S{gIndx, cIndx}.Attenuation;
					% varValue is stimulus type dependent
					switch(class(obj.Stimuli.S{gIndx, cIndx}))
						case 'DW.Wav'
							[t1, t2, t3] = fileparts(obj.Stimuli.S{gIndx, cIndx}.Filename);
							StimInfo.VarValues{n, g} = [t2 t3];
						case 'DW.Noise'
							tmp1 = obj.Stimuli.S{gIndx, cIndx}.LowerFreq;
							tmp2 = obj.Stimuli.S{gIndx, cIndx}.UpperFreq;
							StimInfo.VarValues{n, g} = sprintf('Noise %.1f - %.1f', tmp1, tmp2);
							clear tmp1 tmp2;				
						case 'DW.Tone'
							StimInfo.VarValues{n, g} = obj.Stimuli.S{gIndx, cIndx}.Frequency;
						case 'double'
							StimInfo.VarValues{n, g} = [];
						otherwise 
							warning('%s: unknown stim class %s', mfilename, ...
																		class(obj.Stimuli.S{gIndx, cIndx}));
							StimInfo.VarValues{n, g} = 'unknown';
					end
				end
			end
		end	% END getStimInfo METHOD
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Plotting methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function varargout = plotUnitWaveforms(obj, varargin)
		%------------------------------------------------------------------------
		% Plots overlaid waveforms of spikes
		%------------------------------------------------------------------------
		% H = Data.plotUnitWaveforms
		%	with no arguments, method will plot data from all probes and units
		%	in individual figures
		%	
		% H = Data.plotUnitWaveforms('probe', probenum)
		%	'probe' option will select probenum for plotting
		%
		% H = Data.plotUnitWaveforms('unit', unitnumber)
		%	'unit' selects unit to display. unitnumber must match a unit ID
		%	for the Probes
		%------------------------------------------------------------------------
		
			%------------------------------------------------
			% ensure that probenum is in bounds
			%------------------------------------------------
			% defaults
			probenum = [];
			unitnum = [];
			if ~isempty(varargin)
				a = 1;
				while a <= length(varargin)
					switch upper(varargin{a})	
						case 'PROBE'
							probenum = varargin{a+1};
							if ~between(probenum, 1, length(obj.Probes))
								error('%s: probe must be in range [1:%d]', ...
													mfilename, length(obj.Probes));
							end
							a = a + 2;
						case 'UNIT'
							unitnum = varargin{a+1};
							a = a + 2;
						otherwise
							error('%s: unknown option %s', mfilename, varargin{a});
					end	% END switch
				end	% END while
			end	% END if ~isempty(varargin)
	
			% if probenum is empty, use all probes
			if isempty(probenum)
				probenum = 1:obj.Nprobes;
			end

			% loop through probes
			fIndex = 1;
			for pindx = 1:length(probenum)
				p = probenum(pindx);
				% check unit
				if isempty(unitnum)
					unitnum = obj.Probes(p).cluster;
				end
				for uindx = 1:length(unitnum)
					u = unitnum(uindx);
					% check for matches
					if any(u == obj.Probes(p).cluster)
						H(fIndex) = figure; %#ok<AGROW>
						W = obj.Probes(p).getWaveformsForCluster(u);
						if ~isempty(W)
							[r, c] = size(W);
							tvec = (1000/obj.Probes(p).samprate) * (0:(r-1));
							plot(tvec, W);
							titlestr = { ...
								sprintf('Probe %d Unit %d', p, u), ...
								sprintf('%s', obj.Probes(p).name), ...
								sprintf('%d waveforms', length(obj.Probes(p).t{uindx})), ...
							};
							title(titlestr, 'Interpreter', 'none');
							xlabel('milliseconds');
						end
						fIndex = fIndex + 1;
					else
						fprintf('Unit %d not found for Probe %d\n', u, p);
					end
				end	% END unit loop
			end	% END probe loop
			if nargout
				varargout{1} = H;
			end
		end	% END plotUnitWaveforms
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function varargout = plotRasterAndPSTH(obj, varargin)
		%------------------------------------------------------------------------
		% Plots rasters of grouped data
		%------------------------------------------------------------------------
		% H = Data.plotRasterAndPSTH
		%	with no arguments, method will plot probe 1 with default binwidth
		%	
		% H = Data.plotRasterAndPSTH('probe', probenum)
		%	'probe' option will select probenum for plotting
		%
		% H = Data.plotRasterAndPSTH('binwidth', psth_bin_width)
		%	'binwidth' sets PSTH binsize (in milliseconds)
		%
		% H = Data.plotRasterAndPSTH('unit', unitnumber)
		%	'unit' selects unit to display
		%
		% H = Data.plotRasterAndPSTH('offset', [pretime posttime])
		%	If spikes are desired from before the stim onset timestamp or after 
		%	the next stim onset timestamp, an 'offset' option of form
		%	[pretime posttime] may be included.  times must be in milliseconds
		%------------------------------------------------------------------------
			
			%------------------------------------------------
			% ensure that probenum is in bounds
			%------------------------------------------------
			probenum = 1;
			unitnum = 0;
			plotopts.binwidth = 5;
			offset = [0 0];
			% !!!! this really needs to be determined by the 
			% timestamps, but there's not an easy way at the moment... use
			% this for now......
			plotopts.timelimits = [0 1000];
			if ~isempty(varargin)
				a = 1;
				while a <= length(varargin)
					switch upper(varargin{a})	
						case 'PROBE'
							probenum = varargin{a+1};
							a = a + 2;
						case 'BINWIDTH'
							plotopts.binwidth = varargin{a+1};
							a = a + 2;
						case 'UNIT'
							unitnum = varargin{a+1};
							a = a + 2;
						case 'OFFSET'
							offset = varargin{a+1};
							% make sure pretime is correct sign (should be
							% positive)
							offset(1) = abs(offset(1));
							plotopts.timelimits(1) = plotopts.timelimits(1) - offset(1);
							plotopts.timelimits(2) = plotopts.timelimits(2) + offset(2);
							a = a + 2;
						otherwise
							error('%s: unknown option %s', mfilename, varargin{a});
					end
				end				
			end
			if ~between(probenum, 1, length(obj.Probes))
				error('%s: probe must be in range [1:%d]', mfilename, ...
																	length(obj.Probes));
			else
				fprintf('%s: plotting raster/psth for probe %d, unit %d\n', ...
												mfilename, probenum, unitnum);
			end
			%------------------------------------------------
			% get the spikes struct for probe and unit
			%------------------------------------------------
			S = obj.getSpikesForProbe(probenum, 'unit', unitnum, 'offset', offset);
			%------------------------------------------------
			% loop through groups
			%------------------------------------------------
			ngroups = length(S);
			% initialize H to hold figures
			H = cell(ngroups, 1);
			for g = 1:ngroups
				% get Stimulus List indices for this group
				Sindx = obj.Stimuli.GroupList{g};
				nlevels = length(Sindx);
				allspikes = cell(nlevels, 1);
				plotopts.rowlabels = cell(nlevels, 1);
				% loop through stim indices
				for n = 1:nlevels
					% convert spiketimes to milliseconds
					allspikes{n} = cell(length(S(g).spikes{n}), 1);
					for t = 1:length(S(g).spikes{n})
						allspikes{n}{t} = 0.001*S(g).spikes{n}{t};
					end
					plotopts.rowlabels{n} = sprintf('%d', ...
													obj.Stimuli.S{Sindx(n), 2}.Attenuation);
				end
				% get column labels from stimulus properties
				switch class(obj.Stimuli.S{Sindx(n), 2})
					case 'DW.Noise'
						plotopts.columnlabels{1} = ...
									sprintf('%s: BBN %.1f - %.1f Hz', ...
										obj.fname, ...
										obj.Stimuli.S{Sindx(n), 2}.LowerFreq, ...
										obj.Stimuli.S{Sindx(n), 2}.UpperFreq	);
					case 'DW.Wav'
						plotopts.columnlabels{1} = sprintf('%s: %s', ...
										obj.fname, obj.Stimuli.S{Sindx(n), 2}.Filename);
					case 'DW.Tone'
						plotopts.columnlabels{1} = sprintf('%s: Tone %.1f Hz', ...
										obj.fname, obj.Stimuli.S{Sindx(n), 2}.Freq);
					otherwise
						plotopts.columnlabels{1} = sprintf('%s: %s', ...
										obj.fname, class(obj.Stimuli.S{Sindx(n), 2}));
				end
				figure
				H{g} = rasterpsthmatrix(allspikes, plotopts);
			end	% END g
			if nargout
				varargout{1} = H;
			end
			if nargout == 2
				varargout{2} = S;
			end
		end	% END plotRastersAndPSTH
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function H = plotRasters(obj, varargin)
		%------------------------------------------------------------------------
		% H = Data.plotRasters(probenum)
		%------------------------------------------------------------------------
		% Plots rasters of grouped data.  if probenum is not specified, 
		% Data.Probes(1) will be used!
		%------------------------------------------------------------------------
			
			%------------------------------------------------
			% ensure that probenum is in bounds
			%------------------------------------------------
			probenum = 1;
			unitnum = 0;
			if ~isempty(varargin)
				probenum = varargin{1};
				if ~between(probenum, 1, length(obj.Probes))
					error('%s: probe must be in range [1:%d]', mfilename, ...
																		length(obj.Probes));
				end
				if length(varargin) == 2
					unitnum = varargin{2};
				end
			end
			%------------------------------------------------
			% get the spikes struct for probe 1
			%------------------------------------------------
			S = obj.getSpikesForProbe(probenum, 'unit', unitnum, 'offset', [100 0]);
			%------------------------------------------------
			% loop through groups
			%------------------------------------------------
			ngroups = length(S);
			% initialize H to hold figures
			H = zeros(ngroups, 1);
			for g = 1:ngroups
				H(g) = figure(g);
				% get Stimulus List indices for this group
				Sindx = obj.Stimuli.GroupList{g};
				% loop through stim indices
				for n = 1:length(Sindx)
					subplot(length(Sindx), 1, n)
					% convert spiketimes to milliseconds
					spikes = cell(length(S(g).spikes{n}), 1);
					for t = 1:length(S(g).spikes{n})
						spikes{t} = 0.001*S(g).spikes{n}{t};
					end
					rasterplot(spikes, [-100 1200])
					s = Sindx(n);
					if isa(class(obj.Stimuli.S{s, 2}), 'DW.Wav')
						tstr1 = fullfile(obj.Stimuli.S{s, 2}.Filepath, ...
															obj.Stimuli.S{s, 2}.Filename);
						tstr2 = sprintf('Atten = %d', obj.Stimuli.S{s, 2}.Attenuation);
						title([tstr2 '   ' tstr1], 'Interpreter', 'none');
					end
				end
				xlabel('ms')
			end	% END g		
		end	% END plotRasters
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Spike access methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
			
		%------------------------------------------------------------------------
		function S = getSpikesForProbe(obj, probenum, varargin)
		%------------------------------------------------------------------------
		% S = Data.getSpikesForProbe(probenum)
		% S = Data.getSpikesForProbe(probenum, 'unit', unitnum)
		% S = Data.getSpikesForProbe(probenum, 'offset', [pretime posttime])
		%------------------------------------------------------------------------
		% Stimuli are split by stimulus characteristics (e.g., wav filename,
		% frequency) and then grouped into common values with different 
		% attenuation settings.
		% Spiketimes for all groups in Data.Stimuli.GroupList are retrieved 
		% using this method for one probe
		%
		% If spikes are desired from before the stim onset timestamp or after 
		% the next stim onset timestamp, a second input argument of form
		% [pretime posttime] may be included.  times must be in milliseconds.
		% Note that times will be referenced to the offset start time
		% (stim_onset_timestamp - pretime)
		%
		%------------------------------------------------------------------------
		% S is a (ngroups X 1) struct array with fields:
		%	name			DataWave name for this probe
		%	spikes		[natten x 1] cell array of spiketime data, where each 
		%	 				element is a cell array of [ntrials x 1] vectors of spiketimes
		% 
		% 		to obtain the spikes for sweep (trial) 20 of attenuation  1
		% 		and stimulus group 2, call S(2).spikes{1}{20} :
		% 
		% 			>> S(2).spikes{1}{20}
		% 					ans =
		% 						1.0e+05 *
		% 
		% 						 2.1662
		% 						 2.5055
		% 						 8.2955
		% 
		%		more generally: for sweep SWEEP, atten ATT and stim group STIM:
		%		
		%			S(STIM).spikes{ATT}{SWEEP}
		%------------------------------------------------------------------------
		
			%--------------------------------------------------
			% ensure that probenum is in bounds
			%--------------------------------------------------
			if ~between(probenum, 1, length(obj.Probes))
				error('%s: probe must be in range [1:%d]', mfilename, ...
																	length(obj.Probes));
			end
			%--------------------------------------------------	
			% check varargin
			%--------------------------------------------------
			% default unitnum is 0
			unitnum = 0;
			% set offset to empty
			offset = [];
			% check optional input args
			argn = 1;
			while argn <= length(varargin)
				switch upper(varargin{argn})
					case 'UNIT'
						unitnum = varargin{argn + 1};
						argn = argn + 2;
					case 'OFFSET'
						offset = varargin{argn + 1};
						argn = argn + 2;
					otherwise
						warning('%s: Unknown option %s', mfilename, varargin{argn});
						argn = argn + 2;
				end
			end
			%--------------------------------------------------
			% # of groups available
			%--------------------------------------------------
			ngroups = length(obj.Stimuli.GroupList);
			%--------------------------------------------------
			% get spikes for unit
			%--------------------------------------------------
			S = repmat( struct('spikes', {}, 'name', []), ngroups, 1);
			if isempty(offset)
				% no offset provided
				for g = 1:ngroups
					S(g).spikes = obj.getSpikesForStim(g, probenum, unitnum);
					S(g).name = obj.Probes(probenum).name;
				end
			else
				% use pre/post offset to sweep timestamps
				for g = 1:ngroups
					S(g).spikes = obj.getSpikesForStim(g, probenum, unitnum, ...
																		'offset', offset);
					S(g).name = obj.Probes(probenum).name;
				end
			end	% END if
		end	% END getSpikesForProbe
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function spikes = getSpikesForStim(obj, group, probe, unit, varargin)
		%------------------------------------------------------------------------
		% spikes = Data.getSpikesForStim(group, probe, unit, <<options>>)
		%------------------------------------------------------------------------
		% Stimuli are split by stimulus characteristics (e.g., wav filename,
		% frequency) and then grouped into common values with different 
		% attenuation settings.
		%
		% Spiketimes for each group in Data.Stimuli.GroupList are retrieved 
		% using this method.
		%
		% By default, only spikes between the each Stimulus output timestamp and
		% the following stimulus's output timestamp will be collected.  
		% To extend this time, you may provide an additional 'offset' option
		% is a 1 X 2 vector of pre-stimulus time and post-sweep time
		% in milliseconds:
		%
 		%	spikes = Data.getSpikesForStim(group, probe, unit, 'offset', [pre post])
		%	
		%	e.g.,
		%	
		%	spikes = Data.getSpikesForStim(1, 2, 255, 'offset', [100 200])
		%
		%		This will get spikes for group 1, probe 2, unit id 255 and
		%		include spikes 100 ms before the stimulus sweep timestamp and 
		%		200 ms after the following stimulus timestamp.  If the stimulus
		%		timestamps are, for example, 1000 ms apart, this will include 
		%		spikes in a 1300 ms total window.  
		%
		%	Note that pre-stim onset time is POSITIVE, such that the offset 
		%	sweep timestamp will be computed as:
		%			start_time = original_start_time - pre
		%			end_time = original_end_time + post
		%	or, for our example, 
		%			start_time = original_start_time - 100
		%			end_time = original_end_time + 200
		%
		% Alternatively, you may specify a 'window', relative to the stimulus
		% onset timestamp:
		%
		%	spikes = Data.getSpikesForStim(1, 2, 255, 'window', [-100 575]
		%	
		%	will return spikes that occurred between 100 ms before the stimulus
		%	timestamp and 575 ms after.
		%------------------------------------------------------------------------
		
			%--------------------------------------------------
			% check to make sure group is within bounds
			%--------------------------------------------------
			% # of groups available
			ngroups = length(obj.Stimuli.GroupList);
			if (group < 1) || (group > ngroups)
				error('%s: group must be in range [1:%d]', mfilename, ngroups);
			end
			%--------------------------------------------------
			% check probe number
			%--------------------------------------------------
			if ~between(probe, 1, length(obj.Probes))
				error('%s: probe must be in range [1:%d]', mfilename, ...
																		length(obj.Probes));
			end
			%--------------------------------------------------
			% check unit
			%--------------------------------------------------
			unitchk = (unit == obj.Probes(probe).cluster);
			if ~any(unitchk)
				fprintf('%s: units in probe %d are:\n', mfilename, probe)
				fprintf('\t%d\n', obj.Probes(probe).cluster)
				error('%s: unit %d not found!', mfilename, unit)
			else
				unitid = find(unitchk);
				unitid = unitid(1);
			end
			%--------------------------------------------------
			% check for pre/post time offset or window
			%--------------------------------------------------
			wMode = [];
			if ~isempty(varargin)
				argN = 1;
				while argN <= length(varargin)
					switch upper(varargin{argN})
						case 'OFFSET'
							% assign to pre and post sweeptime (convert to usec)
							presweeptime = 1000 * varargin{argN+1}(1);
							postsweeptime = 1000 * varargin{argN+1}(2);
							wMode = 'OFFSET';
							argN = argN + 2;
						case 'WINDOW'
							% assign to pre and post sweeptime (convert to usec)
							presweeptime = 1000 * varargin{argN+1}(1);
							postsweeptime = 1000 * varargin{argN+1}(2);
							wMode = 'WINDOW';
							argN = argN + 2;
						otherwise
							fprintf('%s.getSpikesForStim: unknown option %s\n', ...
												mfilename, varargin{argN});
							error(mfilename);
					end
				end	% END while argN
			end
			%--------------------------------------------------
			% get the list of stimuli for this group
			%--------------------------------------------------
			Sindx = obj.Stimuli.GroupList{group};
			%--------------------------------------------------
			% figure out window
			%--------------------------------------------------
			% preallocate tstart and tend
			tstart = cell(length(Sindx), 1);
			tend = cell(length(Sindx), 1);
			if isempty(wMode)
				% loop through the groups
				for sloop = 1:length(Sindx)
					s = Sindx(sloop);
					tstart{sloop} = obj.Stimuli.Sweepstart{s};
					tend{sloop} = obj.Stimuli.Sweepend{s};
				end
			elseif strcmpi(wMode, 'OFFSET')
				% loop through the groups and subtract presweeptime from the 
				% sweepstart values and add postsweeptime to the sweepend
				% values
				for sloop = 1:length(Sindx)
					s = Sindx(sloop);
					tstart{sloop} = obj.Stimuli.Sweepstart{s} - abs(presweeptime);
					tend{sloop} = obj.Stimuli.Sweepend{s} + postsweeptime;
				end								
			elseif strcmpi(wMode, 'WINDOW')
				% loop through the groups and subtract presweeptime from the 
				% sweepstart values and add postsweeptime to the sweepstart
				% values
				for sloop = 1:length(Sindx)
					s = Sindx(sloop);
					tstart{sloop} = obj.Stimuli.Sweepstart{s} - abs(presweeptime);
					tend{sloop} = obj.Stimuli.Sweepstart{s} + postsweeptime;
				end
			end	% END if
			%--------------------------------------------------
			% allocate spikes vector
			%--------------------------------------------------
			spikes = cell(length(Sindx), 1);
			%--------------------------------------------------
			% get timestamps (should be in microseconds!)
			%--------------------------------------------------
			Spiketimes = obj.Probes(probe).t{unitid};
			%--------------------------------------------------
			% loop through the groups
			%--------------------------------------------------
			for sloop = 1:length(Sindx)
				s = Sindx(sloop);
				% get the spikes for this stimulus, centered re: Sweepstart
				% timestamp
				spikes{sloop} = find_valid_timestamps(	Spiketimes, ...
																	tstart{sloop}, ...
																	tend{sloop}, ...
																	obj.Stimuli.Sweepstart{s});
			end	% END sloop
			
		end	% END getSpikes
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function out = listProbeInfo(obj)
		%------------------------------------------------------------------------
			if isempty(obj.Nprobes) || ~obj.Nprobes
				warning('%s: Probe data are not available!', mfilename)
				out = [];
				return
			end
			out.name =  cell(obj.Nprobes, 1);
			out.cluster = zeros(obj.Nprobes, 1);
			out.ntimestamps = zeros(obj.Nprobes, 1);

			for n = 1:obj.Nprobes
				out.name{n} = obj.Probes(n).name;
				out.cluster(n) = obj.Probes(n).cluster;
				out.ntimestamps(n) = length(obj.Probes(n).t);
			end
		end	% END listProbeInfo
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
	end	% End of methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
end	% End of classdef
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************

