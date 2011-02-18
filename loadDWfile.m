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

%-----------------------------------------------------------
% get Marker information and retrieve markers from data
%-----------------------------------------------------------
% check for markers
if ~dwinfo.NMarkerCols
	% if no marker cols, error!
	error('%s: no markers detected in file', mfilename)
end

% Pull in Marker data
disp('Reading Probe Data...')

mStartCol = dwinfo.MarkerCols(1);
mEndCol = mStartCol + dwinfo.NMarkerCols;
markerCount = 0;
for L = 1:dwinfo.Ndatalines
	if isempty(data{L}{mStartCol})
		break;
	else
		markerCount = markerCount + 1;

		% save time stamp as number
		Marker(markerCount).t = str2double(data{L}{mStartCol});
		% save times in stand-alone vector
		MarkerTimes = Marker(markerCount).t;
		
		% initialize fields
		Marker(markerCount).text = [];
		Marker(markerCount).value = [];

		% loop through the marker cols
		for M = (mStartCol+1):mEndCol
			tmp = sscanf(data{L}{M}, '%f');
			% check if column datum is a string
			if isempty(tmp)
				% column holds string data
				Marker(markerCount).text = data{L}{M};
			else
				% column hold numeric data
				Marker(markerCount).value = tmp;
			end
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

disp('Reading Probe Data...')
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
