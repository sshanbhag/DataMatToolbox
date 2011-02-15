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

fname = 'Spreadsheet_11Feb11_Format.txt';
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
elseif ~dwinfo.ndata1
	error('%s: no fields found in data line 1 of data file', mfilename);
else
	%-----------------------------------------------------------
	% allocate data cell array
	%-----------------------------------------------------------
	data = cell(dwinfo.Nlines - 2, dwinfo.ndata1);

	%-----------------------------------------------------------
	% open file for text reading
	%-----------------------------------------------------------
	fp = fopen(dwinfo.filename, 'rt');

	%-----------------------------------------------------------
	% loop through first 2 header lines
	%-----------------------------------------------------------
	for n = 1:N_HEADER_LINES
		fgetl(fp);
	end
	
	%-----------------------------------------------------------
	% now, read in data
	%-----------------------------------------------------------
	% set data line counter to 1
	dcount = 1;
	
	% loop through data lines, starting line after header lines
	% (first data line)
	for n = (N_HEADER_LINES + 1):dwinfo.Nlines
		% read in text line from file
		line_n = fgetl(fp);
		% scan in fields 
		tmp = textscan(line_n, '%s', dwinfo.ndata1, 'Delimiter', '\t');
		% check if the number of data cols read in matches with expected #
		ntmp = length(tmp{1});
% 		if ntmp ~= dwinfo.ndata1
% 				warning('%s: data number mismatch in file line %d', mfilename, n)
% 		end
		% loop through the line fields and pull out data strings, store in data
		% cell array for later processing
		for m = 1:ntmp
			data{dcount, m} = tmp{1}{m};
		end
		% increment line counter
		dcount = dcount+1;
	end
	
	%-----------------------------------------------------------
	% close file
	%-----------------------------------------------------------
	fclose(fp);
end

save test.mat data dcount dwinfo

% need to post-process the data a bit to place into suitable vectors




