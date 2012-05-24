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

	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define public properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = public, GetAccess = public)
		file
		path
		Nlines
		header
		firstline = struct('data', {}, 'ndata', []);
		UnitTimestampCols
		MarkerTimestampCols
		MarkerTags
		MarkerCols
		Ncols
		NMarkerCols
		SpikeCols
		NSpikeCols
		UnitCols
		NUnitCols
		Ndatalines
	end
	
	%------------------------------------------------------------------------
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

			obj.file = '';
			obj.path = '';
			LOAD_FLAG = 0;

			%parse input and verify
			if nargin == 0
				return
			else

				errFlag = 0;
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
					end
					n = n + 1;
				end

			end

			%--------------------------------------------------------------
			% Load data
			%--------------------------------------------------------------
			if LOAD_FLAG
				obj.readDataWaveTextInfo;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.file(obj, val)
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
		function [obj, errFlg] = readHeader(obj)
			errFlg = 0;
			
			% first check to make sure filename is set
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
			fprintf('... found %d lines (including header) in file:\n %s \n', obj.Nlines, filename);

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
			
			% parse header if no error
			if ~errFlg
				[~, tmp] = obj.parseHeader;
				if tmp
					errFlg = 3;
				end
			end

		end
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
		%		NMarkerCols			# of marker columns
		%		SpikeCols			indices of spike timestamp columns
		%		NSpikeCols			# of spike columns
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
			%-----------------------------------------------------------
			% load defaults
			%-----------------------------------------------------------
			DataWaveDefaults;

			%--------------------------------------------------------------------
			% find spike timestamp header tags
			%--------------------------------------------------------------------
			% algorithm: search text fields in header using strncmp, then locate
			%  non-zero elements in returned vector using find
			%--------------------------------------------------------------------
			timestamp_cols = strncmp(obj.header.fields{1}, 'timestamp', length('timestamp'));
			obj.SpikeCols = find(timestamp_cols);
			% make sure something was found
			if isempty(obj.SpikeCols)
				% if empty, return error
				errFlg = 2;
				return
			end
			obj.NSpikeCols = length(obj.SpikeCols);

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
			%--------------------------------------------------------------------
			tmp = strncmp(obj.header.fields{1}, 'Timestamp', length('Timestamp'));
			obj.MarkerCols = find(tmp);
			% make sure something was found
			if isempty(obj.MarkerCols)
				% if empty, warn user
				warning('DWFILE:MARKER', '%s: no Marker timestamp fields found in file %s header', ...
												mfilename, obj.filename);
			elseif length(obj.MarkerCols) > 1
				% if unpredicted length, warn user
				warning('DWFILE:MARKER', '%s: %d Marker timestamp  fields found in file %s header', ...
												mfilename, length(obj.MarkerCols), obj.filename);
			else
				% 3 of columns for marker information is the value of the first
				% spike timestamp column - 1 (all spike timestamps are written after
				% the marker information)
				obj.NMarkerCols = obj.SpikeCols(1) - 1;
			end

			%-----------------------------------------------------------
			% get marker tags
			%-----------------------------------------------------------
			obj.MarkerTags = obj.header.fields{1}(1:obj.NMarkerCols);

			obj.Ndatalines = obj.Nlines - N_HEADER_LINES;


		end	%parseHeader
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

	% End of methods
	end
end