function varargout = loadDWfile(varargin)
%------------------------------------------------------------------------
% [D, errFlg, data] = loadDWfile(fname, pname)
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
%	data	"raw" data cell array
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
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% default # of header lines
%-----------------------------------------------------------
N_HEADER_LINES = 2;
N_CHANNELS_PER_PROBE = 5;

%-----------------------------------------------------------
% check input arguments, act depending on inputs
%-----------------------------------------------------------
if nargin == 1
	% only filename was provided, use current directory as path
	fname = varargin{1};
	pname = pwd;
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
[dwinfo, errFlg] = readDWfileinfo(fname, pname);
% check if errFlg ~= 0 (error detected)
if errFlg
	warning('DWFILE:ReadError', '%s readDWfileinfo returned error flag %d', ...
												mfilename, errFlg);
end

%-----------------------------------------------------------
% perform some checks on file information, if okay, read in
% data from DW text file
%-----------------------------------------------------------
if ~dwinfo.header.nfields(2)
	error('%s: no fields found in header line 2 of data file', mfilename);
elseif ~dwinfo.Ncols
	error('%s: no fields found in data line 1 of data file', mfilename);
end

%-----------------------------------------------------------
% allocate data cell array
%-----------------------------------------------------------
%data = cell(dwinfo.Nlines - 2, dwinfo.ndata1);
data = cell(dwinfo.Nlines - N_HEADER_LINES, 1);

%-----------------------------------------------------------
% open file for text reading
%-----------------------------------------------------------
fp = fopen(dwinfo.filename, 'rt');

%-----------------------------------------------------------
% skip past first 2 header lines
%-----------------------------------------------------------
for n = 1:N_HEADER_LINES
	fgetl(fp);
end

%-----------------------------------------------------------
% now, read in data
%-----------------------------------------------------------
disp(['Reading Data from ' dwinfo.filename ' ... ']);
% loop through data lines, starting line after header lines
% (first data line)
for line_index = 1:(dwinfo.Nlines - N_HEADER_LINES)
	% read in text line from file
	line_in = fgetl(fp);
	% scan in fields 
	tmp = textscan(line_in, '%s', dwinfo.Ncols, 'Delimiter', '\t');
	% save in data cell array
	data{line_index} = tmp{1};
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
disp('Parsing Marker Data...')

mStartCol = dwinfo.MarkerCols(1);
mEndCol = mStartCol + dwinfo.NMarkerCols;
markerCount = 0;
for L = 1:dwinfo.Ndatalines
	% check if marker information is present
	if isempty(data{L}{mStartCol})
		% if not, break
		disp(['found end of Marker Data, line ' num2str(L)])
		break;
	else
		% store marker information
		
		% increment marker index
		markerCount = markerCount + 1;

		% save time stamp as number
		Marker(markerCount).t = str2double(data{L}{mStartCol});
		% save times in stand-alone vector
		MarkerTimes(markerCount) = Marker(markerCount).t;
		
		% initialize fields for marker data
		% I am assuming two potential types of data
		%	text		some sort of string, e.g. .wav file name
		%	value		a numeric value (e.g., attenuation value)
		Marker(markerCount).text = [];
		Marker(markerCount).value = [];

		% loop through the marker cols
		m = 1;
		for M = (mStartCol+1):mEndCol
			tmp = sscanf(data{L}{M}, '%f');
			% check if column datum is a string
			if isempty(tmp)
				% column holds string data
				Marker(markerCount).text = data{L}{M};
				Marker(markerCount).vars{m} = data{L}{M};
			else
				% column hold numeric data
				Marker(markerCount).value = tmp;
				Marker(markerCount).vars{m} = tmp;
			end
			m = m+1;
		end
	end
end

%-----------------------------------------------------------
% Pull in Probe Data
%-----------------------------------------------------------

% check # of probes
NumberOfProbes = length(dwinfo.ProbeCols);
if ~NumberOfProbes
	% if 0, error
	error('%s: no probes detected in header', mfilename)
else
	% otherwise, build Probe data structure
	tmpProbe = struct('C0', [], 'C1', [], 'C2', [],'C3', [], 'C4', []);
	for n = 1:NumberOfProbes
		Probe(n) = tmpProbe;
	end
end

disp('Parsing Probe Data...')
% loop through data lines
for l = 1:dwinfo.Ndatalines
	% loop through probes (i.e., tetrodes)
	for p = 1:NumberOfProbes
		% loop through data columns for this probe
		for c = 0:(N_CHANNELS_PER_PROBE-1)
			% column in data cell array 
			datacol = c + dwinfo.ProbeCols(p);
			% if data element at line l, column datacol is not empty,
			% assign the value (convert to double from string) to the
			% appropriate channel in the probe struct array
			if ~isempty(data{l}{datacol})
				Probe(p).(['C' num2str(c)])(l) = str2double(data{l}{datacol});
			end
		end % end c
	end % end p
end % end l

D.info = dwinfo;
D.Probe = Probe;
D.Marker = Marker;
D.MarkerTimes = MarkerTimes;
D.NumberOfProbes = NumberOfProbes;

if any(nargout == [0 1 2 3])
	varargout{1} = D;
end

if any(nargout == [2 3])
	varargout{2} = errFlg;
end

if nargout == 3
	varargout{3} = data;
end
