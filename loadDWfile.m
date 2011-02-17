% function [D, errFlg] = loadDWfile(fname, path)
%------------------------------------------------------------------------
% [D, errFlg] = loadDWfile(fname, path)
%------------------------------------------------------------------------
% 
% Description
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Input		input info
% 
% Output Arguments:
% 	Output	output info
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 24 January, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% default # of header lines
%-----------------------------------------------------------
N_HEADER_LINES = 2;
N_CHANNELS_PER_PROBE = 5;

fname = 'Spreadsheet_17Feb11_Format.txt';
%-----------------------------------------------------------
% get file info, read header
%-----------------------------------------------------------
[dwinfo, errFlg] = readDWfileinfo(fname);
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

save test.mat data n dwinfo
