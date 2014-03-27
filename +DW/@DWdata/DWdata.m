%-----------------------------------------------------------------------------
% DWdata.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW Package
% Class Definition
%-----------------------------------------------------------------------------
%	DWdata facilitates access to data output by loadDWfile.m
%
%  obj = DWdata('xxx_y_zz_ttt_n.mat') loads the following:
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
%	 -	renamed file and class from LoadDWfileData.m to DWdata.m
%	 -	beginning massive Objectification...
%	13 Dec 2012 (SJS):
%	 -	fixing issues on initialization of object
%	14 Jan 2013 (SJS)
%	 -	major overhaul to use NeuroShare Matlab API and DataWave Neuroshare
%		interface object (@NS)
%	31 Jan 2013 (SJS): reworking to avoid segmentation faults 
%	26 Mar 2014 (SJS): seg fault is endemic to Matlab Neuroshare API
%		reworking to use output files from PyDDF project
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
% classdef (ConstructOnLoad = true) DWdata < handle
classdef DWdata < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		fpath			%path to file
		fname			%name of file
		fext			%file extension
		DDF			% neuroshare interface object
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
		function obj = DWdata(varargin)
		%---------------------------------------------------------------------	
		%	DWdata(<fileName>) 
		%	Constructor method
		%	opens file called fileName (char) or opens
		%  a dialog box to get a filename if the fileName provided does not exist.
		%---------------------------------------------------------------------
			% set default input data mode (ddf file)
			mode = 'DDF';
			
			% first, parse input and verify
			if nargin > 1
				% mode input was provided
				fprintf('Using txt file input (from PyDDF)');
				mode = 'TXT';
				
			elseif nargin == 1
				% filename was given, make sure it exists, otherwise throw
				% error
				if exist(varargin{1}, 'file') == 2
					[obj.fpath, obj.fname, obj.fext] = fileparts(varargin{1});
					if strcmpi(obj.fext, 'TXT')
						mode = 'TXT';
					else
						mode = 'DDF';
					end
				else
					error('DWdata: %s does not exist!',varargin{1});
				end
				
			elseif nargin == 0
				% no filename given, open gui panel to get file
				[fileName, obj.fpath] = uigetfile({'*.ddf'; '*.txt'; '*.*'}, ...
																'Open .ddf file from DataWave');
				if fileName == 0
					% if user hit cancel (filename == 0), return
					fprintf('no filename chosen, exiting\n');
					return
				end
				[~, obj.fname, obj.fext] = fileparts(fileName);
				if strcmpi(obj.fext, 'TXT')
						mode = 'TXT';
				else
						mode = 'DDF';
				end
			end
			
			% store in obj.fullfname
			obj.fullfname = fullfile(obj.fpath, [obj.fname obj.fext]);
			
			fprintf('Initializing from file %s...\n\n', obj.fullfname);
			% initialize NS object
			obj.initNS(mode);
		end	% END DWdata CONSTRUCTOR
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Data input/reading methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = initNS(obj, varargin)
		%------------------------------------------------------------------------
		% initNS
		%------------------------------------------------------------------------
		% so long as DW.fullfname is set, initNS method will initialize the
		% DataWave NeuroShare interface object DW.DDF
		%------------------------------------------------------------------------
			if length(varargin)
				mode = varargin{1};
			end
			
			if ~isempty(obj.fullfname) && exist(obj.fullfname, 'file')
				if strcmpi(mode, 'DDF')
					obj.DDF = DW.NS(obj.fullfname);
				else
					obj.DDF = DW.PyDDF(obj.fullfname);
				end
			else
				fprintf('file not found\n')
				return
			end
		end	% END initNS
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function varargout = loadMarkers(obj)
		%------------------------------------------------------------------------
		% [Events, errFlg] = DWdata.loadMarkers
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% Initial setup
			%-----------------------------------------------------------
			DataWaveDefaults	% load defaults
			errFlg = 0;

			%-----------------------------------------------------------
			% read Events using NeuroShare
			%-----------------------------------------------------------
			Events = obj.DDF.getEvents;
			if isempty(Events)
				errFlg = 1;
				warning('%s: error loading Events', mfilename);
				return
			end

			%-----------------------------------------------------------
			% parse Markers
			%-----------------------------------------------------------			
			% create a list of Event IDS for each event entity.
			%	EventID == 0 --> non-marker entity
			%	EventID == 1 --> L channel markers
			%	EventID == 2 --> R channel markers
			EventID = zeros(obj.DDF.nEvent, 1);

			% first, find the R and L marker events (R is usually Events(2), 
			% L is usually Events(3), but try not to assume)
			for v = 1:obj.DDF.nEvent
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
				error('%s: EventCount mismatch!', mfilename);
			else
				EventCount = Events(MarkerEvents(1)).EventCount;
				obj.Nmarkers = EventCount;
			end

			evL = find(EventID == L);
			evR = find(EventID == R);
			fprintf('Left Markers -> Event(%d)\n', evL);
			fprintf('Right Markers -> Event(%d)\n', evR);

			% allocate Markers object array
			% Matlab thinks obj.Markers is a double, and will throw an
			% an error.  so, we trick Matlab in order to do so:  first, need to
			% initialize Markers as a single marker, then allocate 
			% full array.  otherwise, Matlab things
			obj.Markers = DW.Marker;
			obj.Markers(EventCount, 1) = DW.Marker;

			% loop through the events
			for n = 1:EventCount
				% parse into cell arrays, to preserve empty fields, do NOT
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
		% Stimuli = DWdata.loadStimuli
		%------------------------------------------------------------------------
			DataWaveDefaults	% load defaults
			%-----------------------------------------------------------
			% load markers if they haven't been loaded yet
			%-----------------------------------------------------------
			if isempty(obj.Markers)
				obj.loadMarkers
			end
			if ~obj.Nmarkers
				error('%s: no markers in DWdata!')
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
		function loadProbes(obj)
		%------------------------------------------------------------------------
		% Probes = DWdata.loadProbes
		%------------------------------------------------------------------------
			DataWaveDefaults	% load defaults
			%-----------------------------------------------------------
			% load markers if they haven't been loaded yet
			%-----------------------------------------------------------
			if isempty(obj.Markers)
				obj.loadMarkers
			end
			if ~obj.Nmarkers
				error('%s: no markers in DWdata!')
			end
			
			%-----------------------------------------------------------
			% you'd assume that the neural events from 
			% DataWave/Neuroshare are what we want for the spike data.
			% you'd be wrong.  the segment entities have the spiketimes
			% as well as the spike snippets (if stored).
			%-----------------------------------------------------------
			% check if there are any segment entities in the data
			if ~obj.DDF.nSegment
				warning('%s: there are no Segment entities in data!', mfilename);
				return
			end
			
			% load segment timestamps
			Segment = obj.DDF.getSegmentTimeStampsOnly;
	
			obj.Probes = DW.Probe;
			if obj.DDF.nSegment > 1
				obj.Probes(obj.DDF.nSegment, 1) = DW.Probe;
			end
			
			% unit 0 is uncatergorized, 255 is noise
			for n = 1:obj.DDF.nSegment
				 obj.Probes(n).t = Segment(n).TimeStamp(Segment(n).UnitID == 0);
				 obj.Probes(n).cluster = Segment(n).SourceInfo.ProbeInfo;
			end
			clear Segment
		end	% END loadProbes
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

