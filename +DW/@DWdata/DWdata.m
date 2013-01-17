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
%-----------------------------------------------------------------------------
% TO DO: lots....
%-----------------------------------------------------------------------------

classdef (ConstructOnLoad = true) DWdata < handle
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
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
	% Define other properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties
		fullfname;  %full file name
	end
  
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
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

			% first, parse input and verify
			if nargin > 1
				% too many inputs
				error('DWdata: too many inputs!');
				
			elseif nargin == 1
				% filename was given, make sure it exists, otherwise throw
				% error
				if exist(varargin{1}, 'file') == 2
					[obj.fpath, obj.fname, obj.fext] = fileparts(varargin{1});
				else
					error('DWdata: %s does not exist!',varargin{1});
				end
				
			elseif nargin == 0
				% no filename given, open gui panel to get file
				[fileName, obj.fpath] = uigetfile('*.ddf', ...
																'Open .ddf file from DataWave');
				if fileName == 0
					% if user hit cancel (filename == 0), return
					fprintf('no filename chosen, exiting\n');
					return
				end
				[~, obj.fname, obj.fext] = fileparts(fileName);
			end
			
			% store in obj.fullfname
			obj.fullfname = fullfile(obj.fpath, [obj.fname obj.fext]);
			
			fprintf('Initializing from file %s...\n\n', obj.fullfname);
			% initialize NS object
			obj.initNS;
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
		function obj = initNS(obj)
		%------------------------------------------------------------------------
		% initNS
		%------------------------------------------------------------------------
		% so long as DW.fullfname is set, initNS method will initialize the
		% DataWave NeuroShare interface object DW.DDF
		%------------------------------------------------------------------------
			if ~isempty(obj.fullfname) && exist(obj.fullfname, 'file')
				obj.DDF = DW.NS(obj.fullfname);
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

			% read Events using NeuroShare
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
				ts = [Events(evR).TimeStamp(n) Events(evL).TimeStamp(n)];
				if all(ts > 0)
					% check if they are equal
					if ts(1) == ts(2)
						% if so, use the right one
						obj.Markers(n).Timestamp = ts(1);
					else
						% otherwise, store the first one (in time)
						obj.Markers(n).Timestamp = min(ts);
					end
				else
					% store timestamp that is > 0
					obj.Markers(n).Timestamp = ts(ts>0);
				end
				clear ts elist
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
		%------------------------------------------------------------------------
		function varargout = loadStimuli(obj)

			% load markers if they haven't been loaded yet
			if isempty(obj.Markers)
				obj.loadMarkers
			end
			
			if nargout == 0
				return
			else
				varargout{1} = 1;
			end
			
		end	% END loadStimuli
		
		
		
		
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [obj, rawdata, errFlg] = readRawData(obj)
		%------------------------------------------------------------------------
		% [rawdata, rawdata, errFlg] = readRawData
		%------------------------------------------------------------------------
		% formerly readDataWaveTextFile
		% Reads raw text information from Datawave Text file
		% 
		% returns raw text in a 2-D cell array.
		%
		% it is a big hog of memory to store this in the
		% DWdata object itself
		%------------------------------------------------------------------------
		% Output Arguments:
		%
		% 	errFlg	Error flag
		% 					0		no error
		%
		%------------------------------------------------------------------------
		% See: readDataWaveHeader, parseDataWaveTextHeader
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% initial things
			%-----------------------------------------------------------
			errFlg = 0;
			% load defaults
			DataWaveDefaults;

			% check input arguments, act depending on inputs
			if nargin == 0
				% no filename or path provided
				error('%s: input argument error', mfilename);
			end

			% allocate data cell array  
			% * N_HEADER_LINES defined in DataWaveDefaults.m file
			rawdata = cell(obj.Info.Nlines - N_HEADER_LINES, 1);

			%-----------------------------------------------------------
			% Read data file
			%-----------------------------------------------------------
			% open file for text reading
			fp = fopen(obj.fullfname, 'rt');
			% skip past header lines
			for n = 1:N_HEADER_LINES
				fgetl(fp);
			end
			disp(['Reading Data from ' obj.fullfname ' ... ']);
			% loop through rawdata lines, starting line after header lines
			% (first rawdata line)
			for line_index = 1:(obj.Info.Nlines - N_HEADER_LINES)
				% read in text line from file
				line_in = fgetl(fp);
				% scan in fields 
				tmp = textscan(line_in, '%s', obj.Info.Ncols, 'Delimiter', '\t');
				% save in rawdata cell array
				rawdata{line_index} = tmp{1};
			end
			% close file
			fclose(fp);
			% no error(s) encountered
			errFlg = 0;
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function errFlg = loadData(obj)
		%------------------------------------------------------------------------
		% Loads data from file specified in obj.fname.  Initializes
		% Markers, Probes
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% Initial setup
			%-----------------------------------------------------------
			errFlg = 0;
			DataWaveDefaults;
			%-----------------------------------------------------------
			% read raw data
			%	this is done within a method to avoid having to store
			%	the raw DataWave text data in memory or in a file.
			%-----------------------------------------------------------
			[~, rawdata, errFlg] = obj.readRawData;
			if errFlg
				warning('%s: error %d', mfilename, errFlg);
				return
			end
			%-----------------------------------------------------------
			% parse Markers
			%-----------------------------------------------------------
			obj.parseMarkersFromData(rawdata);
			obj.Nmarkers = obj.Info.Nmarkers;
			%-----------------------------------------------------------
			% parse Probes
			%-----------------------------------------------------------
			% determine # of probes from the Info object, NSpikeTimeCols
			% this *should* be accurate, so long as the # of header
			% columns labeled "timestamp" in the DataWave text output
			% file is consistent with the # of probes.
			obj.Nprobes = obj.Info.NSpikeTimeCols;
			% parse the probe information from the raw text data
			errFlg = obj.parseProbesFromData(rawdata);
			if errFlg
				fprintf('loadData: error in DW.parseProbesFromData');
			end
			%-----------------------------------------------------------
			% parse Units
			%-----------------------------------------------------------
			obj.parseProbeIntoUnits;
			%-----------------------------------------------------------
			% clean up
			%-----------------------------------------------------------
			clear rawdata

		end	%END loadData
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function errFlg = parseMarkersFromData(obj, rawdata)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
			%-----------------------------------------------------------
			% Initial setup
			%-----------------------------------------------------------
			errFlg = 0;
			DataWaveDefaults;
			%-----------------------------------------------------------
			% check inputs
			%-----------------------------------------------------------
			if isempty(rawdata)
				error('%s: rawdata is empty!', mfilename)
			elseif ~iscell(rawdata)
				error('%s: improper format for rawdata; must be cell array or vector', mfilename);
			end
			%-----------------------------------------------------------
			% build list of Marker objects
			%-----------------------------------------------------------
			disp( 'Building Marker list ...');
			% get # of markers (# of rawdata values)
			fprintf('%d markers in data\n', obj.Info.Nmarkers)
			% intialize Marker
			obj.Markers = DW.Marker;
			% loop, using # of data lines as upper bound
			for L = 1:obj.Info.Nmarkers
				obj.Markers(L) = DW.Marker(rawdata{L});
				% assign marker #
				obj.Markers(L).N = L;
			end
		end	% END parseMarkersFromData()
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function errFlg = parseProbesFromData(obj, rawdata)
		%------------------------------------------------------------------------
		% errFlg = parseProbesFromData(obj, rawdata)
		%------------------------------------------------------------------------
		
			%-----------------------------------------------------------
			% Initial setup
			%-----------------------------------------------------------
			errFlg = 0;
			DataWaveDefaults;
			
			%-----------------------------------------------------------
			% check inputs
			%-----------------------------------------------------------
			if isempty(rawdata)
				error('%s: rawdata is empty!', mfilename)
			elseif ~iscell(rawdata)
				error('%s: improper format for rawdata; must be cell array or vector', mfilename);
			end		
		
			%-----------------------------------------------------------
			% Pull in Spike Channel Data
			%-----------------------------------------------------------
			obj.Nprobes = obj.Info.NSpikeTimeCols;
	
			% check # of Spike channels
			if ~obj.Nprobes
				% if 0, error
				error('%s: no spike data channels detected in header', mfilename)
			else
				% otherwise, build Probes
				obj.Probes = DW.Probe;
				for n = 1:obj.Nprobes
					obj.Probes(n) = DW.Probe;
				end
			end

			disp('Reading and Parsing Probe Data...')
			% loop through spike columns (i.e., tetrodes)
			for p = 1:obj.Info.NSpikeTimeCols
				% get current column number for spike data
				c = obj.Info.SpikeTimeCols(p);
				% pull out the data for this column
				l = 0;
				tmpline = rawdata{1};
				loopFlag = 1;
				while ~isempty(tmpline{c}) && (l < obj.Info.Ndatalines) && loopFlag
					l = l+1;
					tmpline = rawdata{l};
					if ~isempty(tmpline{c})
						obj.Probes(p).t(l) = str2double(tmpline{c});
						obj.Probes(p).cluster(l) = str2double(tmpline{c+1});
					else
						loopFlag = 0;
					end
				end
			end
		end	% END parseProbesFromData
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function parseProbeIntoUnits(obj)
		%---------------------------------------------------------------------
		
			%-----------------------------------------------------------
			% load defaults
			%-----------------------------------------------------------
			DataWaveDefaults;
			errFlg = 0;

			%-----------------------------------------------------------
			% separate spike times into units
			%-----------------------------------------------------------
			% initialize total unit counter
			NUnits = 0;
			% initialize Units object array
			obj.Units = DW.Unit;
			
			% loop through probes
			for p = 1:obj.Nprobes
				% find unique cluster id values
				unique_clusters = unique(obj.Probes(p).cluster);

				% check if unique values were found
				if ~isempty(unique_clusters)
					% if so, figure out how many clusters there are
					obj.Probes(p).Nclusters = length(unique_clusters);

					% then, assign data to the Units() struct array that will hold
					% unit and timestamp data

					% initialize cluster index n
					n = 1;
					% loop through running count of total units, and store probe, unit number
					% indices and timestamps for this unit
					for u = (1 + NUnits) : (NUnits + obj.Probes(p).Nclusters)
						obj.Units(u) = DW.Unit;
						obj.Units(u).probe = p;
						obj.Units(u).unit = unique_clusters(n);
						if obj.Units(u).unit > 0
							obj.Units(u).sorted = 1;
						else
							obj.Units(u).sorted = 0;
						end
						obj.Units(u).indices = find(obj.Probes(p).cluster == obj.Units(u).unit);
						obj.Units(u).timestamp = obj.Probes(p).t(obj.Units(u).indices);
						n = n + 1;
					end

					NUnits = NUnits + obj.Probes(p).Nclusters;
				else
					obj.Probes(p).Nclusters = 0;
					error('%s: no units found for probe %d', mfilename, p);
				end
			end
			
			obj.Nunits = NUnits;

		end
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------


		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function errFlg = exportRawData(obj, varargin)
		%------------------------------------------------------------------------
		% errFlg = exportRawData(obj, varargin)
		%------------------------------------------------------------------------
		% exportRawData						with no arguments, will ask for a -mat
		%											file to use for output
		%
		% exportRawData(<filename>)
		%											writes DataWave Raw data to file 
		%											specified in <filename> as a -mat file
		%											(default)
		%
		% exportRawData(<filename>, <mode>)
		%		<mode> options:
		%			'text'			writes to <filename> as text (tab-delimited)
		%			'csv'				writes to <filename> as text (comma-delimited)
		%			'mat'				writes to <filename> as mat file
		%------------------------------------------------------------------------
			errFlg = 0;
			
			%-----------------------------------------------------------
			% parse inputs
			%-----------------------------------------------------------
			if isempty(varargin)
				[filen, pathn] = uiputfile(	...
													{	'*.mat', 'MAT-file (*.mat)'; ...
														'*.txt', 'tab-delimited text file (*.txt)'; ...
														'*.csv', 'comma-delimited text file (*.csv)' }, ...
														'select file for output' ...
												);
														
				if isequal(filen, 0) || isequal(pathn, 0)
					fprintf('\n\n%s: export cancelled\n\n', mfilename)
					return
				else
					fullname = fullfile(pathn, filen);
					[~, ~, fmode] = fileparts(fullname);
				end
				
			else
				% user provided at least a filename
				fullname = varargin{1};
				
				% check if other args present
				if length(varargin) > 1
					switch upper(varargin{2})
						case 'TEXT'
							fmode = '.txt';
						case 'CSV'
							fmode = '.csv';
						case 'MAT'
							fmode = '.mat';
						otherwise
							warning('%s: unknown export mode %s', varargin{2});
							errFlg = 1;
							return
					end
				end
			end

			%-----------------------------------------------------------
			% read raw data
			%	this is done within a method to avoid having to store
			%	the raw DataWave text data in memory or in a file.
			%-----------------------------------------------------------
			[~, rawdata, errFlg] = obj.readRawData;
			if errFlg
				warning('%s: error %d', mfilename, errFlag);
				return
			end
			% determine size of rawdata
			[nrows, ~] = size(rawdata);
			ncols = length(rawdata{1});
			
			%-----------------------------------------------------------
			% write raw data
			%-----------------------------------------------------------
			switch fmode
				% write text file
				case {'.txt', '.csv'}
					% set delimiter
					if strcmp(fmode, '.txt')
						dlm = sprintf('\t');
					elseif strcmp(fmode, '.csv')
						dlm = ',';
					end
					% open file
					fp = fopen(fullname, 'w');
					% write rawdata
					for r = 1:nrows
						for c = 1:ncols
							fprintf(fp, '%s%c', rawdata{r}{c}, dlm);
						end
						fprintf(fp, '\n');
					end
					% close file
					fclose(fp);
					
				case '.mat'
					% save to MAT file
					save(fullname, 'rawdata', '-MAT');
					
				otherwise
					% uh-oh
					error('%s: unknown fmode %s', mfilename, fmode);
			end
			%-----------------------------------------------------------
			% clean up
			%-----------------------------------------------------------
			clear rawdata
		
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------	
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function arrOut = getUniqueWavNames(obj, channel)
		%------------------------------------------------------------------------
		% arrOut = getUniqueWavNames(obj, channel)
		%------------------------------------------------------------------------
		% checks for different .wav file stimulus files
		%	Nunique			# of unique strings in Marker.text
		%	uniqueNames		unique values of strings in Marker.text, cell array
		%	uniqueIndices	indices where each of the unique strings are found in
		%						Marker.WavFilename(L/R)
		%
		% if list of names is empty, leave wavFiles struct empty
		%------------------------------------------------------------------------
			if isnumeric(channel)
				if channel == 1
					fieldname = 'WavFilenameL';
				else
					fieldname = 'WavFilenameR';
				end
			else
				fieldname = ['WavFilename' upper(channel)];
			end
			
			allwavnames = obj.getMarkerFieldAsArray(fieldname);
			
			[names, indices, N] = findUniqueText(allwavnames);
			if ~isempty(names)
				arrOut = struct(	'uniqueNames', names, ...
										'uniqueIndices', indices, ...
										'Nunique', N ...
									);
			else
				arrOut = [];
			end
			
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function arrOut = getMarkerFieldAsArray(obj, fieldname)
		%------------------------------------------------------------------------
		% arrOut = getMarkerFieldAsArray(obj, fieldname)
		% DWdata class method
		%------------------------------------------------------------------------
		% for a given <fieldname> that is a valid marker field  
		% (see DataWaveDefaults for values of MARKER_TAGS), returns an 
		% [Nmarkers X 1] cell vector of values	(if field refers to a string or
		% vector) or a [Nmarkers X 1] vector of values (numerical)
		%------------------------------------------------------------------------
			if ~exist('fieldname', 'var')
				error('%s: no fieldname provided', mfilename)
			end
			% load defaults
			DataWaveDefaults;			
			% check if fieldname is valid
			fcompare = strcmp(fieldname, MARKER_TAGS); %#ok<USENS>
			if isempty(fcompare)
				error('%s: %s is not a valid marker tag field', mfilename, fieldname);
			end
			
			% find index location of fieldname in MARKER_TAGS for dynamic indexing
			floc = find(fcompare);
			% return arrOut as cell array for character fields
			if strcmpi(MARKER_TYPES{floc}, 'char') %#ok<USENS>
				arrOut = cell(obj.Info.Nmarkers, 1);
				for c = 1:obj.Info.Nmarkers
					arrOut{c} = obj.Markers(c).(MARKER_TAGS{floc});
				end
		
			else
				% return as numerical array for number fields
				arrOut = zeros(obj.Info.Nmarkers, 1);
				for c = 1:obj.Info.Nmarkers
					arrOut(c) = obj.Markers(c).(MARKER_TAGS{floc});
				end
			end
		end	% END getMarkerFieldAsArray
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%{
		% Set/Get Methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.fname(obj, val)
			% check if no filename or path provided
			if isempty(val)
				obj.fname = '';
				return
			elseif isnumeric(val)
				warning('%s: fname must be a string; file property unchanged', mfilename);
				return
			else
				% assume path is included or that file is in current directory
				[obj.fpath, obj.fname, obj.fext] = fileparts(val);
				obj.fullfname = fullfile(obj.fpath, [obj.fname obj.fext]);
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.fpath(obj, val)
			% check if no filename or path provided
			if isempty(val)			
				obj.fpath = '';
				return
			elseif isnumeric(val)
				warning('%s: fpath must be a string; path property unchanged', mfilename);
				return
			else
				% user provided path string
				obj.fpath = val;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function ret = get.fullfname(obj)
			ret = fullfile(obj.fpath, [obj.fname obj.fext]);
		end	%get.fullfname
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		%}
		
	% End of methods
	end
end