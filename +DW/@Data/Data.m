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
% 		DDF			% neuroshare interface object
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

			% first, parse input and verify
			if nargin == 0
				return
			end
			% D struct was provided
			obj.initFromStruct(varargin{1});
			% check if a filename was provided
			if nargin > 1
				obj.fullfname = varargin{2};
				[obj.fpath, obj.fname, obj.fext] = fileparts(varargin{2});
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
		
			% store info bits in obj.Info!
			infotags = {'FileType', 'EntityCount', 'TimeStampResolution', ...
								'TimeSpan', 'AppName', 'Time_Year', 'Time_Month', ...
								'Time_Day', 'Time_Hour', 'Time_Min', 'Time_Sec', ...
								'Time_MilliSec', 'FileComment' };
			for n = 1:length(infotags)
				obj.Info.(infotags{n}) = D.(infotags{n});
			end
% 			% store other bits in obj.DDF!
% 			ddftags = {'EntityInfo', 'EventList', 'AnalogList', 'SegmentList', ...
% 							'NeuralList', 'nNeural', 'nSegment', 'nAnalog', 'nEvent'};
% 			for n = 1:length(ddftags)
% 				obj.DDF.(ddftags{n}) = D.(ddftags{n});
% 			end
			% store Markers
			obj.loadMarkers(D.Event);
			% convert to stimuli
			obj.loadStimuli;
			% convert segments to Probes
			obj.loadProbes(D.Segment);
			% get sweep times
			obj.Stimuli.extractTimeFromMarkers(obj.Markers);
			% find stimulus groups (same apart from attenuation)
			obj.Stimuli.findCommon;
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
			% is stored in an individual Event entity: the left (channel
			% 1) marker is in EventID = 1, the right (channel 2) is in 
			% EventId 2.
			% 
			% create a list of Event IDS for each event entity.
			% 
			%	EventID == 0 --> non-marker entity
			%	EventID == 1 --> L channel markers
			%	EventID == 2 --> R channel markers
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

			% check to ensure that both L and R Entities were found
			if ~all( any(EventID == L) & any(EventID == R) )
				% if not abort
				error('%s: Markers for both (L and R) channels not found!', mfilename);
			end
			% check MarkerEvents
			MarkerEvents = find(EventID);
			nMarkerEvents = length(MarkerEvents);
			if nMarkerEvents ~= 2
				error('%s: nMarkerEvents is %d (must be equal to 2)',...
								mfilename, nMarkerEvents);
			end
			% make sure # of events are the same
			if ~(Events(MarkerEvents(1)).EventCount == Events(MarkerEvents(2)).EventCount)
				error('%s: EventCount mismatch between L and R channels!', mfilename);
			else
				EventCount = Events(MarkerEvents(1)).EventCount;
				obj.Nmarkers = EventCount;
			end

			evL = find(EventID == L);
			evR = find(EventID == R);
			fprintf('Left Markers -> Event(%d)\n', evL);
			fprintf('Right Markers -> Event(%d)\n', evR);

			%-----------------------------------------------------------			
			% allocate Markers object array
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
		function varargout = loadProbes(obj, Segment)
		%------------------------------------------------------------------------
		% Probes = Data.loadProbes
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
				nSegments = length(Segment);
			end
		
			%-----------------------------------------------------------
			% you'd assume that the neural events from 
			% DataWave/Neuroshare are what we want for the spike data.
			% you'd be wrong.  the segment entities have the spiketimes
			% as well as the spike snippets (if stored).
			%-----------------------------------------------------------
			% create the Probe object
			obj.Probes = DW.Probe;
			if nSegments > 1
				obj.Probes(nSegments, 1) = DW.Probe;
			end	
			
			% pull out timestamps based on unit ID
			% according to Neuroshare, unit 0 is uncatergorized, 255 is noise
			unitID = cell(nSegments, 1);
			for n = 1:nSegments
				% find unique unit IDs
				unitID{n} = unique(Segment(n).UnitID);
				% store timestamps for units with ID 0 - should
				% probably figure out a way to store/sort all unitIDs at some
				% point........
				obj.Probes(n).t = Segment(n).TimeStamp(Segment(n).UnitID == 0);
				obj.Probes(n).cluster = Segment(n).SourceInfo.ProbeInfo;
			end
			clear Segment
			obj.Nprobes = nSegments;
			if nargout
				varargout{1} = unitID;
			end
		end	% END loadProbes
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
		function spikes = getSpikesForStimulusGroup(obj, groupnum, probenum)
			% # of groups available
			ngroups = length(obj.Stimuli.GroupList);
			% check to make sure groupnum is within bounds
			if ~between(groupnum, 1, ngroups)
				error('%s: group must be in range [1:%d]', mfilename, ngroups)
			end
			if ~between(probenum, 1, length(obj.Probes))
				error('%s: probe must be in range [1:%d]', mfilename, ...
																		length(obj.Probes));
			end
			% get the list of stimuli for this group
			Sindx = obj.Stimuli.GroupList{groupnum};
			% allocate spikes vector
			spikes = cell(length(Sindx), 1);
			% get timestamps (in seconds!)
			Spiketimes = obj.Probes(probenum).t;
			% convert to microseconds
			Spiketimes = 1e6 * Spiketimes;
			% loop through the groups
			for sloop = 1:length(Sindx)
				s = Sindx(sloop);
				% get the spikes for this stimulus
				spikes{sloop} = find_valid_timestamps(Spiketimes, ...
																  obj.Stimuli.Sweepstart{s}, ...
																	obj.Stimuli.Sweepend{s});
			end
		end	% END getSpikes
		
		
		
		
	end	% End of methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
end	% End of classdef
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************

