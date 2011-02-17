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
N_CHANNELS_PER_PROBE = 5;

%-----------------------------------------------------------
% load test.mat file containing variables:	data{}  dcount dwinfo
%-----------------------------------------------------------
if ~exist('data', 'var')
	load test.mat
end

% get number of rows and cols in Data
[Drows, Dcols] = size(data);


%-----------------------------------------------------------
% get Marker information and retrieve markers from data
%-----------------------------------------------------------
% check for markers
if ~dwinfo.NMarkerCols
	% if no marker cols, error!
	error('%s: no markers detected in file', mfilename)
end

% Pull in Marker data
mStartCol = dwinfo.MarkerCols(1);
mEndCol = mStartCol + dwinfo.NMarkerCols;

markerCount = 0;
for L = 1:Drows
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

tic
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
method1 = toc




