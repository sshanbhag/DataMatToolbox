%-----------------------------------------------------------------------------
% DWinfo.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% Class Definition
%-----------------------------------------------------------------------------
% Properties
% 		file						filename
% 		path						path to file = pname;
% 		Nlines					# of lines (including header) in file
% 		header					header structure
%				line					cell vector of header line text
%				fields				cell vector of field names
%				nfields				# of tab-delimited fields in each header line
%		firstline				struct to store first line of data - used to
%									determine fields, spikes, etc.
%				data					text from data line 1
%				ndata					# of fields in data line 1
%		Nmarkers
%		UnitTimestampCols		list of data columns with time stamp data
%		MarkerTimestampCols	list of data columns with marker stime stamp data
% 		MarkerTags				variable names in Marker data
% 		MarkerCols				indices of marker columns
% 		NMarkerCols				# of marker columns
% 		SpikeCols				indices of spike timestamp columns
% 		NSpikeCols				# of spike columns
% 		UnitCols					indices to unit id columns
% 		NUnitCols				# unit id columns (length(UnitCols))
% 		Ndatalines				# of data lines (adjusted for # of header lines)
% 									dwinfo.Nlines - N_HEADER_LINES
%-----------------------------------------------------------------------------
% See also: Marker, DWdata, loadDWfile (function)
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 21 May, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

classdef (ConstructOnLoad = true) DWinfo < handle
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------

	
	%% ----------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define public properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = public, GetAccess = public)
		file
		path
		Nlines
		Ndatalines
		header
		firstline = struct('data', {}, 'ndata', []);
		Nmarkers
		Ncols
		MarkerTags
		MarkerCols
		NMarkerCols
		MarkerTimeCols
		NMarkerTimeCols
		SpikeTimeCols
		NSpikeTimeCols
		UnitCols
		NUnitCols
	end
	
	%% ----------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = DWinfo(varargin)
		%-----------------------------------------------------------------	
		%	DWinfo(<fileName>) 
		%	Constructor method
		%-----------------------------------------------------------------

			%--------------------------------------------------------------
			%	Set file/path information
			%--------------------------------------------------------------

			% determine optional argument length by checking size of varargin.
			Noptarg = size(varargin, 2);
			% determine "std" args by taking difference of total # args and Noptarg
			Nstdarg = nargin - Noptarg;
			
			%--------------------------------------------------------------
			%parse input and verify
			%--------------------------------------------------------------
			obj.file = '';
			obj.path = '';
			LOAD_FLAG = 0;
			errFlag = 0;
			
			% check inputs
			if nargin == 0
				return
				
			else
				n = 1;
				
				% loop through Noptarg
				while (n <= Noptarg) && ~errFlag
					
					% check for file option string
					if strcmpi(varargin{n}, 'File')
						% increment index
						n = n + 1;
						% check if outta bounds
						if n > Noptarg
							error('%s: no argument provided for option %s', mfilename, varargin{n-1});
						end
						% assume path is included or that file is in current directory 
						[pname, fname, ext] = fileparts(varargin{n});
						if isempty(pname)
							% no path, assume current dir.
							obj.path = pwd;
						else
							obj.path = pname;
						end
						obj.file = [fname ext];
						clear pname fname ext
						n = n + 1;

					% check for Load option 
					elseif strcmpi(varargin{n}, 'Load')
						% increment index
						n = n + 1;
						LOAD_FLAG = 1;

					% unknown input
					else
						warning('%s: unknown option %s', mfilename, varargin{n});
						n = n + 1;
					end
					
				end

			end

			%--------------------------------------------------------------
			% Load data
			%--------------------------------------------------------------
			if LOAD_FLAG
				fprintf('%s: reading and parsing Header...\n', mfilename);
				obj.readAndParseHeader;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.file(obj, val)
		%-----------------------------------------------------------------	
		%	file set method
		%-----------------------------------------------------------------

			% check if no filename or path provided
			if isempty(val)			
				obj.file = '';
				return
			elseif isnumeric(val)
				warning('%s: file must be a string; file property unchanged', mfilename);
				return
			else
				% assume path is included or that file is in current directory
				[pname, fname, ext] = fileparts(val);
				% set values for path and file
				obj.path = pname;
				obj.file = [fname ext];
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.path(obj, val)
			% check if no filename or path provided
			if isempty(val)			
				obj.path = '';
				return
			elseif isnumeric(val)
				warning('%s: path must be a string; path property unchanged', mfilename);
				return
			else
				% user provided path string
				obj.path = val;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = openFile(obj, varargin)

			% check if no filename
			if size(varargin, 2) == 0			
				% open ui panel to get filename and path from user
				[fname, pname] = uigetfile('*.txt', 'Select Datawave Exported Text File', obj.path);
				% return if user pressed cancel button
				if isequal(fname, 0) || isequal(pname,0)
					disp('Cancelled...')
					return
				else
					% pull fname and ext from the fname returned by uigetfile
					% pname is already defined from uigetfile, so ignore that here
					[~, fname, ext] = fileparts(fname);
				end
			else 
				% assume path is included or that file is in current directory
				[pname, fname, ext] = fileparts(varargin{1});
			end

			% set the path
			obj.path = pname;
			% save filename
			obj.file = [fname ext];
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = openPath(obj, varargin)

			% check if no filename or path provided
			if size(vararing, 2) == 0			
				% open ui panel to get filename and path from user
				pname = uigetdir('Select Datawave Data Directory', obj.path);
				% return if user pressed cancel button
				if isequal(pname,0)
					disp('Cancelled...')
					return
				end
			else
				% user provided path string
				pname = varargin{1};
			end

			obj.path = pname;
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [obj, errFlg] = readAndParseHeader(obj)
		%------------------------------------------------------------------------
		% [obj, errFlg] = readAndParseHeader(obj)
		%------------------------------------------------------------------------

			% first check to make sure filename is set
			errFlg = 0;
			if isempty(obj.file)
				warning('%s: file is undefined', mfilename);
				errFlg = 1;
				return
			end
		
			%-----------------------------------------------------------
			% read header if no error
			%-----------------------------------------------------------
			[~, errFlg] = obj.readHeader;
			if errFlg
				fprintf('%s: skipping parsing of header due to error %d\n', mfilename, errFlg);
				return
			end
			
			%-----------------------------------------------------------
			% parse header if no error
			%-----------------------------------------------------------
			[~, errFlg] = obj.parseHeader;
			if errFlg
				fprintf('%s: error %d while parsing header\n', mfilename, errFlg);
			end
			
			%-----------------------------------------------------------
			% count markers if no error
			%-----------------------------------------------------------
			[~, nmarkers, errFlg] = obj.countMarkers;
			if errFlg
				fprintf('%s: error %d while counting markers\n', mfilename, errFlg);
			else
				fprintf('...found %d markers in file\n', nmarkers)
			end
			
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [obj, errFlg] = readHeader(obj)
		%------------------------------------------------------------------------
		% [obj, errFlg] = readHeader(obj)
		%------------------------------------------------------------------------

			% first check to make sure filename is set
			errFlg = 0;
			if isempty(obj.file)
				warning('%s: file is undefined', mfilename);
				errFlg = 1;
				return
			end

			% load defaults
			DataWaveDefaults;

			% pre-allocate header struct
			obj.header = struct(	'line',  cell(N_HEADER_LINES, 1), ...
									'fields', cell(N_HEADER_LINES, 1), ...
									'nfields', zeros(N_HEADER_LINES, 1)		);

			% count the number of lines in the file
			filename = fullfile(obj.path, obj.file);
			obj.Nlines = countTextFileLines(filename);
			fprintf('... found %d lines (including header) in file:\n\t\t %s \n', ...
								obj.Nlines, filename);

			% # of data lines are total number of lines in the 
			% file minus the # of header lines specified in DataWaveDefaults file
			obj.Ndatalines = obj.Nlines - N_HEADER_LINES;
			
			% Open file for reading as text
			fp = fopen(filename, 'rt');

			% Read in header lines and parse to get file information
			for n = 1:N_HEADER_LINES
				obj.header.line{n} = fgetl(fp);
				% get values in line n of header (text, tab-delimited)
				tmp = textscan(obj.header.line{n}, '%s', 'Delimiter', '\t');
				obj.header.fields{n} = tmp{1};
				obj.header.nfields(n) = length(obj.header.fields{n});
			end

			% warn if no fields counted
			if ~sum(obj.header.nfields)
				warning('DWFILE:empty-data', ...
							'%s: no fields found in header lines of data file', mfilename);
				errFlg = 2;
			end

			%-----------------------------------------------------------
			% read in 1st data line and count fields
			%-----------------------------------------------------------
			tmpline = fgetl(fp);
			tmp = textscan(tmpline, '%s', 'Delimiter', '\t');
			obj.firstline(1).data = tmp{1};
			obj.firstline.ndata = length(obj.firstline.data);
			obj.Ncols = obj.firstline.ndata;

			%-----------------------------------------------------------
			% close file
			%-----------------------------------------------------------
			fclose(fp);

		end  % end readHeader
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------


		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [obj, errFlg] = parseHeader(obj)
		%------------------------------------------------------------------------
		% [obj, errFlg] = parseHeader
		%------------------------------------------------------------------------
		% 	Parses information from header to find:
		%		MarkerCols			indices of marker columns
		%		NMarkerCols			# of marker colums
		%		MarkerTimeCols		indices of marker timestamp columns
		%									(should be 1 in most cases)
		%		NMarkerTimeCols	# of marker timestamp columns
		%		SpikeTimeCols		indices of spike timestamp columns
		%									(should start at 34 usually)
		%		NSpikeTimeCols		# of spike timestamp columns
		%		UnitCols				indices to unit id columns
		%		NUnitCols			# unit id columns (length(UnitCols))
		%		Ndatalines			# of data lines (adjusted for # of header lines)
		%								dwinfo.Nlines - N_HEADER_LINES
		%		MarkerTags			updated MarkerTags
		% 
		%------------------------------------------------------------------------
		% Input Arguments:
		%	none
		% 
		% Output Arguments:
		%	obj		self
		% 
		% 	errFlg	Error flag
		% 					0		no error
		% 					1		user cancelled file opening
		% 					2		no fields found in header lines
		% 					3		file not found
		%------------------------------------------------------------------------
		% See: readDataWaveTextInfo, loadDWfile
		%------------------------------------------------------------------------
			errFlg = 0;
			if isempty(obj.file)
				warning('%s: file is undefined', mfilename);
				errFlg = 1;
				return
			end
			
			%-----------------------------------------------------------
			% load defaults
			%-----------------------------------------------------------
			DataWaveDefaults;

			%--------------------------------------------------------------------
			% find spike timestamp header tags
			%--------------------------------------------------------------------
			% algorithm: search for 'timestamp' in header using strncmp, 
			% then locate non-zero elements in returned vector using find
			%--------------------------------------------------------------------
			timestamp_cols = strncmp(obj.header.fields{1}, 'timestamp', length('timestamp'));
			obj.SpikeTimeCols = find(timestamp_cols);
			% make sure something was found
			if isempty(obj.SpikeTimeCols)
				% if empty, return error
				errFlg = 2;
				return
			end
			% store # timestamp cols
			obj.NSpikeTimeCols = length(obj.SpikeTimeCols);
			
			%--------------------------------------------------------------------
			% find unit header tags in header
			%--------------------------------------------------------------------
			tmp = strncmp(obj.header.fields{1}, 'cluster', 1);
			obj.UnitCols = find(tmp);
			% make sure something was found
			if isempty(obj.UnitCols)
				% if empty, warn user
				warning('DWFILE:TSTAMP', '%s: no unit/probe timestamp fields found in file %s header', ...
												mfilename, obj.filename);
			end
			obj.NUnitCols = length(obj.UnitCols);

			%--------------------------------------------------------------------
			% find marker Timestamp header tags
			%--------------------------------------------------------------------
			% There are several ways to do this.  
			% one that makes very few assumptions is to AND together the unit and spike
			% timestamp cols, NOT them and those should be the marker cols.
			%			- OR	-
			% search for 'Timestamp'
			%--------------------------------------------------------------------
			tmp = strncmp(obj.header.fields{1}, 'Timestamp', length('Timestamp'));
			obj.MarkerTimeCols = find(tmp);
			% make sure something was found
			if isempty(obj.MarkerTimeCols)
				% if empty, warn user
				warning('DWFILE:MARKER', '%s: no Marker timestamp fields found in file %s header', ...
												mfilename, obj.filename);
			elseif length(obj.MarkerTimeCols) > 1
				% if unpredicted length, warn user
				warning('DWFILE:MARKER', '%s: %d Marker timestamp  fields found in file %s header', ...
												mfilename, length(obj.MarkerTimeCols), obj.filename);
			else
				% 3 of columns for marker information is the value of the first
				% spike timestamp column - 1 (all spike timestamps are written after
				% the marker information)
				obj.NMarkerCols = obj.SpikeTimeCols(1) - 1;
				obj.MarkerCols = 1:obj.NMarkerCols;
				obj.NMarkerTimeCols = length(obj.MarkerTimeCols);
			end

			%-----------------------------------------------------------
			% get marker tags
			%-----------------------------------------------------------
			% marker tags are fields in 1st line up to NMarkerCols
			obj.MarkerTags = obj.header.fields{1}(1:obj.NMarkerCols);
			
		end	%parseHeader
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [obj, markerCount, errFlg] = countMarkers(obj)
		%------------------------------------------------------------------------
		% [obj, markerCount, errFlg] = countMarkers
		%------------------------------------------------------------------------
		% figure out # of markers
		%------------------------------------------------------------------------
		% algorithm:
		% 	read in lines of data after header, keeping running count
		% 	when first empty marker timestamp (first element
		% 	in cell array of read text fields) is found, break and 
		% 	store running count as # of markers
		%------------------------------------------------------------------------

			%-----------------------------------------------------------
			% make sure file exists
			%-----------------------------------------------------------
			errFlg = 0;
			if isempty(obj.file)
				warning('%s: file is undefined', mfilename);
				errFlg = 1;
				return
			end
			
			%-----------------------------------------------------------
			% load defaults
			%-----------------------------------------------------------
			DataWaveDefaults;
			
			%-----------------------------------------------------------
			% now, count lines
			%-----------------------------------------------------------
			% Open file for reading as text
			fp = fopen( fullfile(obj.path, obj.file), 'rt');

			% Marker timestamp column - hardcoded here, but should figure out
			% a way to make this programmatically determined.
			mStartCol = 1;

			% initialize marker Counter variable
			markerCount = 0;

			% skip header lines
			for n = 1:N_HEADER_LINES
				fgetl(fp);			
			end
			
			% loop, using # of data lines as upper bound
			for L = 1:obj.Ndatalines
				% read line
				tmpline = fgetl(fp);
				% get values in line L of header (text, tab-delimited)
				tmpvals = textscan(tmpline, '%s', 'Delimiter', '\t');
				% check if marker information is present
				if isempty(tmpvals{1}{mStartCol})
					% if not, break
					break;
				else
					% store marker information
					% increment marker index
					markerCount = markerCount + 1;
				end
			end
			
			% store marker count in Nmarkers property
			obj.Nmarkers = markerCount;
		
		end % end count markers
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

	% End of methods
	end
end