function [out, errFlg] = readDataWaveTextInfo(varargin)
%------------------------------------------------------------------------
% [out, errFlg] = readDataWaveTextInfo(fname, pname)
%------------------------------------------------------------------------
% 
% Gets information from a DataWave exported text file
% 
% Reads information from the the 1st row of text fields (header line)
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	fname		name of DataWave exported text file (usually with .txt file
% 				extension)
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
% 
% 		data1						text from data line 1
% 		ndata1					# of fields in data line 1
%		UnitTimestampCols		list of data columns with time stamp data
%		MarkerTimestampCols	list of data columns with marker stime stamp data
% 		MarkerTags				variable names in Marker data
%
% 	errFlg	Error flag
% 					0		no error
% 					1		user cancelled file opening
% 					2		no fields found in header lines
% 					3		file not found
%
%------------------------------------------------------------------------
% See: readDataWaveHeader 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 25 January, 2011 (SJS)
%
% Revisions:
%	14 Feb 2011 (SJS):
% 		-	adjustments made for new file format with n probes, no
% 			"spacing" columns between probes and markers
% 	27 Apr 2011 (SJS):
% 		-	renamed TimestampCols to UnitTimestampCols - better reflects
% 			data stored there
% 			created MarkerTimestampCols to indicate column with marker timestamps
%	19 May, 2011 (SJS):
% 		NO LONGER USEFUL WITH DataWave txt files due to format change 
% 	20 May, 2011 (SJS): renamed to readDataWaveTextInfo.m
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
disp(['... found ' num2str(Nlines) ' lines in file ' fname ' (including header).']);

%-----------------------------------------------------------
% Open file for reading as text
%-----------------------------------------------------------
fp = fopen(filename, 'rt');

%-----------------------------------------------------------
% Read in header lines and parse to get file information
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
% assign values to output structure
%-----------------------------------------------------------
out.file = fname;
out.path = pname;
out.filename = filename;
out.Nlines = Nlines;
out.header = header;
out.data1 = data1;
out.Ncols = ndata1;

% no error(s) encountered
errFlg = 0;

