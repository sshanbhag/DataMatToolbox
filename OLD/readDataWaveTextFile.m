function [rawdata, errFlg] = readDataWaveTextFile(dwinfo)
%------------------------------------------------------------------------
% [rawdata, errFlg] = readDataWaveTextFile(dwinfo)
%------------------------------------------------------------------------
% 
% Reads raw text information from Datawave Text file
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	dwinfo	DataWave file information object 
% 				(from name of DataWave exported text file (usually with .txt file
% 				extension)
%				If fname is not provided, a GUI window will open to select file.
% 
% Output Arguments:
% 	rawdata		{Nlines X 1} cell array of string vectors corresponding to 
% 					lines of text from file dwinfo.filename
%
% 	errFlg	Error flag
% 					0		no error
%
%------------------------------------------------------------------------
% See: readDataWaveHeader, parseDataWaveTextHeader
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 14 July, 2011 (SJS)
%
% Revisions:
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
if nargin == 0
	% no filename or path provided
	error('%s: input argument error', mfilename);
end

%-----------------------------------------------------------
% allocate data cell array  
% * N_HEADER_LINES defined in DataWaveDefaults.m file
%-----------------------------------------------------------
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
% read in raw data using textscan - this will load
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

% no error(s) encountered
errFlg = 0;

