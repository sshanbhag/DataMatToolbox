function validtimes = find_valid_timestamps(timestamps, winstart, winend)
%------------------------------------------------------------------------
% validtimes = find_valid_timestamps(timestamps, winstart, winend)
%------------------------------------------------------------------------
% locates timestamps that are within the time window given by 
% winstart and winend.  winstart and winend can be vectors, but must
% be of the same length!
% For proper behaviour, please ensure that time units for timestamps
% and windows are the same.
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
	fprintf('Warning %s: no timestamps in %s\n', mfilename, inputname(1))
	validtimes = {};
	return
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
		validtimes{w} = timestamps(validSpikes) - winstart(w);
	else
		validtimes{w} = [];
	end
end
