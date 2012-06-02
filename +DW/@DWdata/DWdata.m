%-----------------------------------------------------------------------------
% DWdata.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
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
% Continuing development: 
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: January, 2012 (TS)
%
% Revisions:
%	17 May, 2012 
%		updated documentation
%		renamed file and class from LoadDWfileData.m to DWdata.m
%		beginning massive Objectification...
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

classdef (ConstructOnLoad = true) DWdata < handle
	
	%% ------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		fpath;			%path to file
		fname;			%name of file
		fext;			%file extension
		Info;			
		Background;	%really this should be a class too...
		D;				%really this should be a class too...
		Stimulus;		%really this should be a class too...
		Markers
		Probes
		Units
		Nmarkers
		Nprobes
		Nunits
	end
	
	%% ------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define other properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties
		fullfname;  %full file name
	end
  
	%% ------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods
	
		%% Constructor
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = DWdata(varargin)
		%---------------------------------------------------------------------	
		%	DWdata(<fileName>) 
		%	Constructor method
		%	opens file called fileName (char) or opens
		%  a dialog box to get a filename if the fileName provided does not exist.
		%---------------------------------------------------------------------	

			%parse input and verify
			if nargin > 1
				error('DWdata:toomanyinputs','too many inputs! Try DWdata(filename)');
			elseif nargin == 1
				if exist(varargin{1}, 'file') == 2
					[obj.fpath, obj.fname, obj.fext] = fileparts(varargin{1});
				else
					warning('DWdata:sillyfname','%s does not exist!',varargin{1});
					return
				end
			elseif nargin == 0
				return;
			end

			%if still don't have a file name, get one
			if isequal(obj.fname,'')
				[fileName, obj.fpath] = uigetfile('*.txt','Open .txt file outupt from DataWave');
				if fileName == 0
					error('DWdata:nochoose','no filename chosen')
				end
				[~, obj.fname, obj.fext] = fileparts(fileName);
			end

			fprintf('Initializing from file %s...\n\n', obj.fullfname);
			
			obj.Info = DW.DWinfo('File', obj.fullfname, 'Load');
			
			%{
			%Try loading the file and see if crashes
			try
				loaded = load(obj.fullfname, 'D', 'Stimulus', 'Background');
				obj.D = loaded.D;
				obj.Stimulus = loaded.Stimulus;
				obj.Background = loaded.Background;
			catch
				error('%s: could not find file %s', mfilename, obj.fullfname);
			end
			%}
		end		%DWdata
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%% General Methods
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function disp(obj)
		%------------------------------------------------------
		% displays information about loaded data
		%------------------------------------------------------
			fprintf(1, '\t%s - Loaded from file:\n',class(obj));
			fprintf(1,'\t\t%s\n\n',obj.fullfname);
			fprintf(1,'Contains data:\n');
			fprintf('\tStimulus:\n\t\t%ix%i %s\n',size(obj.Stimulus,1),size(obj.Stimulus,2),class(obj.Stimulus))
			fprintf('\tBackground:\n\t\t%ix%i %s\n',size(obj.Background,1),size(obj.Background,2),class(obj.Background))
			fprintf('\tD:\n\t\t%ix%i %s\n',size(obj.D,1),size(obj.D,2),class(obj.D))
		end	%disp
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		
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
		function errFlg = loadMarkers(obj)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% Initial setup
			%-----------------------------------------------------------
			errFlg = 0;
			DataWaveDefaults;
			%-----------------------------------------------------------
			% read raw data
			%-----------------------------------------------------------
			[~, rawdata, errFlag] = obj.readRawData;
			if errFlag
				warning('%s: error %d', mfilename, errFlag);
				return
			end
			%-----------------------------------------------------------
			% parse Markers
			%-----------------------------------------------------------
			obj.parseMarkersFromData(rawdata);
			obj.Nmarkers = obj.Info.Nmarkers;
		end	%END loadMarkers
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function errFlg = loadData(obj)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% Initial setup
			%-----------------------------------------------------------
			errFlg = 0;
			DataWaveDefaults;
			%-----------------------------------------------------------
			% read raw data
			%-----------------------------------------------------------
			[~, rawdata, errFlg] = obj.readRawData;
			if errFlg
				warning('%s: error %d', mfilename, errFlag);
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
			obj.Nprobes = obj.Info.NSpikeTimeCols;
			errFlg = obj.parseProbesFromData(rawdata)
			%-----------------------------------------------------------
			% parse Units
			%-----------------------------------------------------------
			obj.parseProbeIntoUnits;
			
		end	%END loadMarkers
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
				fieldname = ['WavFilename' upper(channel)]
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
		%------------------------------------------------------------------------
			if ~exist('fieldname', 'var')
				error('%s: no fieldname provided', mfilename)
			end
			
			DataWaveDefaults;
			
			% check if fieldname is valid
			fcompare = strcmp(fieldname, MARKER_TAGS);
			
			if isempty(fcompare)
				error('%s: %s is not a valid marker tag field', mfilename, fieldname);
			end
			
			% find index location of fieldname in MARKER_TAGS for dynamic indexing
			floc = find(fcompare);
			
			% return arrOut as cell array for character fields
			if strcmpi(MARKER_TYPES{floc}, 'char')
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
		

		%% Set/Get Methods
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
		
	% End of methods
	end
end