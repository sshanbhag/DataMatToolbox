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
%										line		cell vector of header line text
%										fields	cell vector of field names
%										nfields	# of tab-delimited fields in each 
% 													header line
% 		data1						text from data line 1
% 		ndata1					# of fields in data line 1
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
	properties (SetAccess = private, GetAccess = public)
		file;
		path;
		Nlines;
		header;
		data1;
		ndata1;
		UnitTimestampCols;
		MarkerTimestampCols;
		MarkerTags;
		MarkerCols;
		NMarkerCols;
		SpikeCols;
		NSpikeCols;
		UnitColsL;
		NUnitCols;
		Ndatalines;
	end
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define public properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties
		testprop;
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
			obj.set.path('');
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
	% 	function obj = setFile(obj, varargin)
		function set.file(obj, val)

			% check if no filename or path provided
			if isempty(val)			
				% open ui panel to get filename and path from user
				[fname, pname] = uigetfile('*.txt', 'Select Datawave Exported Text File');
				% return if user pressed cancel button
				if isequal(fname, 0) || isequal(pname,0)
					disp('Cancelled...')
					obj.file = '';
					return
				else
					[~, fname, ext] = fileparts(fname);
				end
			else 
				% assume path is included or that file is in current directory
				[pname, fname, ext] = fileparts(val);
			end

			if isempty(pname)
				% no path, assume current dir.
				obj.set.path(pwd);
			else
				obj.set.path(pname);
			end
			obj.file = [fname ext];

			obj
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.path(obj, val)

			% check if no filename or path provided
			if isempty(val)			
				% open ui panel to get filename and path from user
				pname = uigetdir('Select Datawave Data Directory');
				% return if user pressed cancel button
				if isequal(pname,0)
					disp('Cancelled...')
					return
				end
			else
				% user provided path string
				pname = varargin{1};
			end

			if isempty(pname)
				% no path, assume current dir.
				obj.path = '';
			else
				obj.path = pname;
			end
	% 		obj.set.file([fname ext]);

			obj
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
	% 	function obj = setFile(obj, varargin)
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
			obj.set.path(pname);
			% save filename
			obj.set.file([fname ext]);

			obj
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

			obj.set.path(pname);

			obj
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = readDataWaveTextInfo(obj)

			% first check to make sure filename is set
			if isempty(obj.file)
				error('%s: file is undefined', mfilename)
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
			disp(['... found ' num2str(obj.Nlines) ' lines in file ' filename ' (including header).']);

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
			obj.data1 = tmp{1};
			obj.ndata1 = length(obj.data1);
			obj.Ncols = ndata1;

			%-----------------------------------------------------------
			% close file
			%-----------------------------------------------------------
			fclose(fp);

		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------


		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function parseString(obj)

		end	%parseString
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

	% End of methods
	end
end