function varargout = loadDWfile(varargin)
%------------------------------------------------------------------------
% [D, errFlg, rawdata] = loadDWfile(fname, pname)
%------------------------------------------------------------------------
% 
% Description
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	fname		name of DataWave exported text file (usually with .dat file
% 				extension.
%				If fname is not provided, a GUI window will open to select file.
% 
% 	pname		path to DataWave file
%				If pname is not provided, current directory will be assumed
% 
% Output Arguments:
% 	D			output structure
%
% 		info		DataWave file information structure:
%			file		filename
%			path		path to file = pname;
%			Nlines	# of lines (including header) in file
%			header	header structure
%							line			cell vector of header line text
%							fields		cell vector of field names
%							nfields		# of tab-delimited fields in each header line
%			data1		text from data line 1
%			ndata1	# of fields in data line 1
% 
% 		Marker
% 		
%		MarkerTimes
% 		
% 		Probe
% 
% 	errFlg	Error flag
% 					0		no error
% 					1		user cancelled file opening
% 					2		no fields found in header lines
% 					3		file not found
%
%	rawdata	"raw" data cell array
%------------------------------------------------------------------------
% See also: readDWfileinfo 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 24 January, 2011 (SJS)
%
% Revisions:
%	17 Feb 2011 (SJS)
% 	 - finished parsing of data{} into Probe and Marker structs
% 	 -	converted to function
% 	 -	documentation
%	27 Apr 2011 (SJS): added some user-feedback bits
%	13 May 2011 (SJS): fixed problem with use of fully-spec'ed filenames
%  25 May, 2011 (SJS): adapting for new data format
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%-----------------------------------------------------------
% check input arguments, act depending on inputs
%-----------------------------------------------------------
if nargin == 1
	% only filename was provided, assume path is included or that file is
	% in current directory 
	[pname, fname, ext] = fileparts(varargin{1});
	if isempty(pname)
		pname = pwd;
	end
	fname = [fname ext];
elseif nargin == 2
	% filename and path provided as input
	fname = varargin{1};
	pname = varargin{2};
elseif nargin == 0
	% no filename or path provided
	% open ui panel to get filename and path from user
	[fname, pname] = uigetfile('*.txt', 'Select Datawave Exported Text File');
	% return if user pressed cancel button
	if isequal(fname, 0) || isequal(pname,0)
		disp('Cancelled...')
		for n = 1:nargout
			varargout{n} = [];
		end
		varargout{2} = 1;
		return
	end
else
	error('%s: input argument or file error', mfilename);
end


%-----------------------------------------------------------
% get file info, read header
%-----------------------------------------------------------
[dwinfo, errFlg] = readDataWaveTextInfo(fname, pname);
% check if errFlg ~= 0 (error detected)
if errFlg
	warning('DWFILE:ReadError', '%s readDataWaveTextInfo returned error flag %d', ...
												mfilename, errFlg);
end

%-----------------------------------------------------------
% perform some checks on file information, if okay, read in
% data from DW text file
%-----------------------------------------------------------
if ~dwinfo.header.nfields(1)
	save('loadDWfile.error.mat', 'dwinfo', '-MAT');
	error('%s: no fields found in header line 2 of data file', mfilename);
elseif ~dwinfo.Ncols
	save('loadDWfile.error.mat', 'dwinfo', '-MAT');
	error('%s: no fields found in data line 1 of data file', mfilename);
end

%-----------------------------------------------------------
% parse header data to get information about data set
%-----------------------------------------------------------
[dwinfo, errFlg] = parseDataWaveTextHeader(dwinfo);

%-----------------------------------------------------------
% allocate data cell array
%-----------------------------------------------------------
%data = cell(dwinfo.Nlines - 2, dwinfo.ndata1);
rawdata = cell(dwinfo.Nlines - N_HEADER_LINES, 1);

%-----------------------------------------------------------
% open file for text reading
%-----------------------------------------------------------
fp = fopen(dwinfo.filename, 'rt');

%-----------------------------------------------------------
% skip past header lines
%-----------------------------------------------------------
for n = 1:N_HEADER_LINES
	fgetl(fp);
end

%-----------------------------------------------------------
% now, read in rawdata using textscan - this will load
% the entire file into a cell array
%-----------------------------------------------------------
disp(['Reading Data from ' dwinfo.filename ' ... ']);
% loop through rawdata lines, starting line after header lines
% (first rawdata line)
for line_index = 1:(dwinfo.Nlines - N_HEADER_LINES)
	% read in text line from file
	line_in = fgetl(fp);
	% scan in fields 
	tmp = textscan(line_in, '%s', dwinfo.Ncols, 'Delimiter', '\t');
	% save in rawdata cell array
	rawdata{line_index} = tmp{1};
end

%-----------------------------------------------------------
% close file
%-----------------------------------------------------------
fclose(fp);

%-----------------------------------------------------------
% get Marker information and retrieve markers from data
%-----------------------------------------------------------
% check for markers
if ~dwinfo.NMarkerCols
	% if no marker cols, error!
	error('%s: no markers detected in file', mfilename)
end

% Pull in Marker data
disp('Reading Marker Data...')

mStartCol = dwinfo.MarkerCols(1);
mEndCol = dwinfo.NMarkerCols;
markerCount = 0;
for L = 1:dwinfo.Ndatalines
	% check if marker information is present
	if isempty(rawdata{L}{mStartCol})
		% if not, break
		disp(['found end of Marker Data, line ' num2str(L)])
		break;
	else
		% store marker information
		
		% increment marker index
		markerCount = markerCount + 1;

		% save time stamp as number
		MarkerData(markerCount).t = str2double(rawdata{L}{mStartCol});
		% save times in stand-alone vector
		MarkerTimes(markerCount) = MarkerData(markerCount).t;
		
		% initialize fields for marker data
		% I am assuming two potential types of data
		%	text		some sort of string, e.g. .wav file name
		%	value		a numeric value (e.g., attenuation value)
		MarkerData(markerCount).string = rawdata{L}(mStartCol:mEndCol);
	end
end



%-----------------------------------------------------------
% Pull in Spike Channel Data
%-----------------------------------------------------------
% check # of Spike channels
if ~dwinfo.NSpikeCols
	% if 0, error
	error('%s: no spike data channels detected in header', mfilename)
else
	% otherwise, build Probe data structure
	for n = 1:dwinfo.NSpikeCols
		Probe(n) = struct('t', [], 'cluster', []);
	end
end

disp('Parsing Probe Data...')
% loop through spike columns (i.e., tetrodes)
for p = 1:dwinfo.NSpikeCols
	% get current column number for spike data
	c = dwinfo.SpikeCols(p);
	% pull out the data for this column
	l = 0;
	tmpline = rawdata{1};
	loopFlag = 1;
	while ~isempty(tmpline{c}) && (l < dwinfo.Ndatalines) && loopFlag
		l = l+1;
		tmpline = rawdata{l};
		if ~isempty(tmpline{c})
			Probe(p).t(l) = str2double(tmpline{c});
			Probe(p).cluster(l) = str2double(tmpline{c+1});
		else
			loopFlag = 0;
		end
	end
end




D.info = dwinfo;
D.Probe = Probe;
D.Marker = parseDataWaveMarkers(MarkerData, dwinfo);
D.MarkerTimes = MarkerTimes;

if any(nargout == [0 1 2 3])
	varargout{1} = D;
end

if any(nargout == [2 3])
	varargout{2} = errFlg;
end

if nargout == 3
	varargout{3} = rawdata;
end
