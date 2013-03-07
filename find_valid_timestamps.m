function varargout = find_valid_timestamps(timestamps, winstart, winend, varargin)
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
% By default, the returned timestamps will be referenced to the time 
% provided in winstart.  This is done under the assumption that each
% winstart value is the start of a "sweep" or stimulus presentation and
% spiketimes are best interpreted relative to that timestamp.
%
% To reference times to an alternative "zero" time, provide winref as input
% (in same units as timestamps, please).
%
% ¡The timestamps may be returned unaltered by using 0 as the winref value!
%
%%------------------------------------------------------------------------
% Input Arguments:
% 	timestamps		vector of timestamps
%	winstart			vector (or single value) of window start time
%	winend			vector (or single value) of window end time
%
%	Optional:
% 		winref		vector (or single value) of reference time
% 		
% Output Arguments:
%	validtimes		cell array of valid timestamps.  cell array will have
%						length equal to # of winstart (and winend) values
%
%	validindices	cell array of valid indices
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
%	6 Mar 2013 (SJS): added some info to help docs.
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
	% otherwise, use the provided value as reference time
	winref = varargin{1};
	if length(winref) == 1
		% if winref is a single value, create a vector the length of the
		% number of time windows for use in each "sweep"
		winref = winref * ones(size(winstart));
	end
end

%------------------------------------------------------------------------
% get values
%------------------------------------------------------------------------
% allocate spikes vector
validtimes = cell(nwin, 1);
validindices = validtimes;
% loop through the groups
for w = 1:nwin
	% find spikes that are within the current window
	validSpikes = (winstart(w) < timestamps) & (timestamps < winend(w));
	validindices{w} = validSpikes;
	% store spiketimes (relative to sweep start)
	if ~isempty(validSpikes)
		validtimes{w} = timestamps(validSpikes) - winref(w);
	else
		validtimes{w} = [];
	end
end

%------------------------------------------------------------------------
% assign outputs
%------------------------------------------------------------------------
varargout{1} = validtimes;
if nargout == 2
	varargout{2} = validindices;
end