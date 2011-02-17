function [out, errFlg] = readDWfileinfo(varargin)
%------------------------------------------------------------------------
% [out, errFlg] = readDWfileinfo(fname, pname)
%------------------------------------------------------------------------
% 
% returns structure 
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
% 	out	DataWave file information structure:
% 		file		filename
% 		path		path to file = pname;
% 		Nlines	# of lines (including header) in file
% 		header	header structure
%						line			cell vector of header line text
% 						fields		cell vector of field names
% 						nfields		# of tab-delimited fields in each header line
% 		data1		text from data line 1
% 		ndata1	# of fields in data line 1
% 
% 	errFlg	Error flag
% 					0		no error
% 					1		user cancelled file opening
% 					2		no fields found in header lines
% 					3		file not found
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 25 January, 2011 (SJS)
%
% Revisions:
%	14 Feb 2011 (SJS):
% 		-	adjustments made for new file format with n probes, no
% 			"spacing" columns between probes and markers
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% default # of header lines
%-----------------------------------------------------------
N_HEADER_LINES = 2;

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
		out = [];
		errFlg = 1;
		return
	end
else
	error('%s: input argument or file error', mfilename);
end

%-----------------------------------------------------------
% initialize filename (full path) and errFlg
%-----------------------------------------------------------
filename = fullfile(pname, fname);
if ~exist(filename, 'file')
	errFlg = 3;
	out = [];
	return
else
	errFlg = 0;
end

%-----------------------------------------------------------
% count the number of lines in the file
%-----------------------------------------------------------
Nlines = countTextFileLines(filename);
disp(['Found ' num2str(Nlines) ' lines in file ' fname ' (including header).']);

%-----------------------------------------------------------
% Open file for reading as text
%-----------------------------------------------------------
fp = fopen(filename, 'rt');

%-----------------------------------------------------------
% Read in 2 header lines and parse to get file information
%-----------------------------------------------------------
header.line = cell(N_HEADER_LINES, 1);
header.fields = cell(N_HEADER_LINES, 1);
header.nfields = zeros(N_HEADER_LINES, 1);
for n = 1:N_HEADER_LINES
	header.line{n} = fgetl(fp);
	% get values in line n of header (text, tab-delimited)
	tmp = textscan(header.line{n}, '%s', 'Delimiter', '\t');
	header.fields{n} = tmp{1};
	header.nfields(n) = length(header.fields{n});
end

% warn if no fields counted
if ~sum(header.nfields)
	warning('DWFILE:empty-data', ...
				'%s: no fields found in header lines of data file', mfilename);
	errFlg = 2;
end

%-----------------------------------------------------------
% read in 1st data line and count fields
%-----------------------------------------------------------
tmpline = fgetl(fp);
tmp = textscan(tmpline, '%s', 'Delimiter', '\t');
data1 = tmp{1};
ndata1 = length(data1);

%-----------------------------------------------------------
% close file
%-----------------------------------------------------------
fclose(fp);

%-----------------------------------------------------------
% parse header information
%-----------------------------------------------------------

% find probe header tags
tmp = strncmp(header.fields{1}, 'Probe', 5);
ProbeCols = find(tmp);
% make sure something was found
if isempty(ProbeCols)
	% if empty, return error
	errFlg = 2;
	out = [];
	return
end
Nprobes = length(ProbeCols);

% find Marker header tags
tmp = strncmp(header.fields{1}, 'Marker', 6);
MarkerCols = find(tmp);
% make sure something was found
if isempty(MarkerCols)
	% if empty, warn user
	warning('DWFILE:MARKER', '%s: no Marker fields found in file %s header', ...
									mfilename, filename);
elseif length(MarkerCols) > 1
	% if unpredicted length, warn user
	warning('DWFILE:MARKER', '%s: %d Marker fields found in file %s header', ...
									mfilename, length(MarkerCols), filename);
else
	NMarkerCols = header.nfields(2) - MarkerCols(1);
end

% find timestamp header tags
tmp = strncmp(header.fields{2}, 'timestamp', 1);
TimestampCols = find(tmp);
% make sure something was found
if isempty(TimestampCols)
	% if empty, warn user
	warning('DWFILE:TSTAMP', '%s: no probe timestamp fields found in file %s header', ...
									mfilename, filename);
end

%-----------------------------------------------------------
% assign values to output structure
%-----------------------------------------------------------
out.file = fname;
out.path = pname;
out.filename = filename;
out.Nlines = Nlines;
out.header = header;
out.data1 = data1;
out.Ncols = ndata1;
out.ProbeCols = ProbeCols;
out.MarkerCols = MarkerCols;
out.NMarkerCols = NMarkerCols;
out.TimestampCols = TimestampCols;
out.Nprobes = Nprobes;
out.Ndatalines = Nlines - N_HEADER_LINES;


