% function Output = function_template(Input)
%------------------------------------------------------------------------
% Output = function_template(Input)
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
N_COLS_PER_PROBE = 6;
N_CHAN_PER_PROBE = 5;

%-----------------------------------------------------------
% load test.mat file containing variables:	data{}  dcount dwinfo
%-----------------------------------------------------------
if ~exist('data', 'var')
	load test.mat
end

[Drows, Dcols] = size(data);


%-----------------------------------------------------------
% first, parse # of probes
%-----------------------------------------------------------
% find columns in header line 1 fields that match 'Probe'
% these will indicate the start column for an N-Probe
tmp = strncmpi('Probe', dwinfo.header.fields{1}, 5);
ProbeCols = find(tmp);
NumberOfProbes = length(ProbeCols);

%-----------------------------------------------------------
% then, find the marker col
%-----------------------------------------------------------
tmp = strncmpi('Marker', dwinfo.header.fields{1}, 6);
MarkerCols = find(tmp);

%-----------------------------------------------------------
% find column headers
%-----------------------------------------------------------
% Timestamp columns
tmp = strncmpi('Timestamp', dwinfo.header.fields{2}, 6);
TimestampCols = find(tmp);
% Other columns
OtherCols = find(~tmp);

%-----------------------------------------------------------
% get Marker information and retrieve markers from data
%-----------------------------------------------------------
if length(MarkerCols) == 1
	tmp = dwinfo.header.fields{2}(MarkerCols:end)
else
	error('%s: strange, length(MarkerCols) ~= 1', mfilename);
end
MarkerNames = cell(length(tmp), 1);
for n = 1:length(tmp)
	MarkerNames{n} = tmp(n);
end

% Pull in Marker data
nMarkers = length(MarkerNames);
mStartCol = MarkerCols(1);
markerCount = 0;
for l = 1:Drows
	if isempty(data{l, mStartCol})
		break;
	else
		markerCount = markerCount + 1;
		Marker.t(l) = str2double(data{l, mStartCol});
		for m = 1:nMarkers-1
			Marker.mval{l, m} = data{l, mStartCol + m};
		end
	end
end

%-----------------------------------------------------------
% Pull in Probe Data
%-----------------------------------------------------------
% check # of probes
if NumberOfProbes == 0
	% if 0, error
	error('%s: no probes detected in header', mfilename)
else
	% otherwise, build Probe data structure
	tmpProbe = struct('C0', [], 'C1', [], 'C2', [],'C3', [], 'C4', []);
	for n = 1:NumberOfProbes
		Probe(n) = tmpProbe;
	end
end

Probe2 = Probe;

tic
% loop through data lines
for l = 1:Drows
	% loop through probes (i.e., tetrodes)
	for p = 1:NumberOfProbes
		% loop through data columns for this probe
		for c = 0:(N_CHAN_PER_PROBE-1)
			% column in data cell array 
			datacol = c + ProbeCols(p);
			% if data element at line l, column datacol is not empty,
			% assign the value (convert to double from string) to the
			% appropriate channel in the probe struct array
			if ~isempty(data{l, datacol})
				Probe(p).(['C' num2str(c)])(l) = str2double(data{l, datacol});
			end
		end % end c
	end % end p
end % end l
method1 = toc


tic
% loop through data lines
for l = 1:Drows
	% loop through probes (i.e., tetrodes)
	for p = 1:NumberOfProbes
		% loop through data columns for this probe
		for c = 0:(N_CHAN_PER_PROBE-1)
			% column in data cell array 
			datacol = c + ProbeCols(p);
			
			% if data element at line l, column datacol is not empty,
			% assign the value (convert to double from string) to the
			% appropriate channel in the probe struct array
			if ~isempty(data{l, datacol})
				Probe(p).(['C' num2str(c)])(l) = str2double(data{l, datacol});
			end
		end % end c
	end % end p
end % end l
method2 = toc



