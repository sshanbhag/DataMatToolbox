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

%-----------------------------------------------------------
% load test.mat file containing variables:	data{}  dcount dwinfo
%-----------------------------------------------------------
if ~exist('data', 'var')
	load test.mat
end

[nrows, ncols] = size(data);


%-----------------------------------------------------------
% first, parse # of probes
%-----------------------------------------------------------
% find columns in header line 1 fields that match 'Probe'
% these will indicate the start column for an N-Probe
tmp = strncmpi('Probe', dwinfo.header.fields{1}, 5);
ProbeCols = find(tmp);

%-----------------------------------------------------------
% then, find the marker col
%-----------------------------------------------------------
tmp = strncmpi('Marker', dwinfo.header.fields{1}, 6);
MarkerCols = find(tmp);

%-----------------------------------------------------------
% find column headers
%-----------------------------------------------------------
tmp = strncmpi('Timestamp', dwinfo.header.fields{2}, 6);
TimestampCols = find(tmp);

%-----------------------------------------------------------
% get Marker information
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


%-----------------------------------------------------------
% Pull in Marker info
%-----------------------------------------------------------
nMarkers = length(MarkerNames);
mStartCol = MarkerCols(1);
markerCount = 0;
for l = (N_HEADER_LINES + 1):length(data)
	if isempty(data{l, mStartCol})
		break;
	else
		markerCount = markerCount + 1;
		Marker.t(l) = data{l, mStartCol};
		for m = 1:nMarkers-1
			Marker.mval{l, m} = data{l, mStartCol + m};
		end
	end
	
	
	
end





