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
		Info = DW.DWinfo;			
		Background;	%really this should be a class too...
		D;				%really this should be a class too...
		Stimulus;		%really this should be a class too...
		Markers
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

			obj.Info.file = obj.fullfname;
			obj.Info.readHeader;
		
			
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
		function ret = get.fullfname(obj)
			ret = fullfile(obj.fpath, [obj.fname obj.fext]);
		end	%get.fullfname
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
		function parseMarkers(obj, rawdata)
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
			% get Marker information and retrieve markers from data
			%-----------------------------------------------------------
			disp( 'Processing Marker information ...')

			% get # of markers (# of rawdata values)
			Nmarkers = length(rawdata);
			fprintf('%d markers in data\n', Nmarkers)

			% initialize variables
			markerCount = 0;
			obj.Markers = DW.Marker;
			
			% loop, using # of data lines as upper bound
			for L = 1:Nmarkers
				% check if marker information is present
				if isempty(rawdata{L})
					% if not, break
					disp(['found end of Marker Data, line ' num2str(L)])
					break;
				else
					% increment marker index
					markerCount = markerCount + 1
					% parse and store marker information
					obj.Markers(markerCount) = DW.Marker(rawdata{L});
					% assign marker #
					obj.Markers(markerCount).N = markerCount;
				end
			end
			
			%{
			%-----------------------------------------------------------------------------
			% Now check for different .wav file stimulus files
			%-----------------------------------------------------------------------------
			% Nunique			# of unique strings in Marker.text
			% uniqueNames		unique values of strings in Marker.text, cell array
			% uniqueIndices	indices where each of the unique strings are found in
			%						Marker.WavFilename(L/R)
			%
			% if list of names is empty, leave wavFiles struct empty
			%-----------------------------------------------------------------------------
			% R channel
			[names, indices, N] = findUniqueText(obj.WavFilenameR);
			% check if names is empty
			if ~isempty(names)
				obj.Marker.wavFilesR = struct(	'uniqueNames', names, ...
													'uniqueIndices', indices, ...
													'Nunique', N ...
													);
			end
			% L channel
			[names, indices, N] = findUniqueText(obj.Marker.WavFilenameL);
			if ~isempty(names)
				obj.Marker.wavFilesL = struct(	'uniqueNames', names, ...
													'uniqueIndices', indices, ...
													'Nunique', N ...
													);
			end
			%}			
		
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
	% End of methods
	end
end