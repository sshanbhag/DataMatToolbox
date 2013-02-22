function validtimes = find_valid_timestamps(timestamps, winstart, winend, varargin)
%------------------------------------------------------------------------
% validtimes = find_valid_timestamps(timestamps, winstart, winend)
% validtimes = find_valid_timestamps(timestamps, winstart, winend, winref)
%------------------------------------------------------------------------
% locates timestamps that are within the time window given by 
% winstart and winend.  winstart and winend can be vectors, but must
% be of the same length!  
% For proper behaviour, please ensure that time units for timestamps
% and windows are the same.
%
% To reference times to an alternative "0" time, provide winref as input
% (in same units as timestamps, please).
%
% returns cell array validtimes
%%------------------------------------------------------------------------
% Input Arguments:
% 	timestamps		vector of timestamps
%	winstart			vector (or single value) of window start time
%	winend			vector (or single value) of window end time
%
% Output Arguments:
%	validtimes		cell array of valid timestamps.  cell array will have
%						length equal to # of winstart (and winend) values
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 12 February, 2013 (SJS)
%
% Revisions:
%	21 Feb 2013 (SJS): changed return validtimes to {[]} when 
%								timestamps is empty for consistency
%	22 Feb 2013 (SJS): added optional input winref
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Check inputs
%------------------------------------------------------------------------
% make sure # of winstart values == # winend values
if length(winstart) ~= length(winend)
	error('%s: need same number of window starts as window ends!', mfilename);
else
	nwin = length(winstart);
end
% make sure there are some timestamps to check!
if isempty(timestamps)
% 	fprintf('Warning %s: no timestamps in %s\n', mfilename, inputname(1))
	validtimes = {[]};
	return
end
if isempty(varargin)
	% if no winref provided, use default winstart
	winref = winstart;
else
	winref = varargin{1};
	if length(winref) == 1
		winref = winref * ones(size(winstart));
	end
end

%------------------------------------------------------------------------
% get values
%------------------------------------------------------------------------
% allocate spikes vector
validtimes = cell(nwin, 1);
% loop through the groups
for w = 1:nwin
	% find spikes that are within the current window
	validSpikes = (winstart(w) < timestamps) & (timestamps < winend(w));
	% store spiketimes (relative to sweep start)
	if ~isempty(validSpikes)
		validtimes{w} = timestamps(validSpikes) - winref(w);
	else
		validtimes{w} = [];
	end
end
