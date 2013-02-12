function validtimes = find_valid_timestamps(timestamps, winstart, winend)

if length(winstart) ~= length(winend)
	error('%s: need same number of window starts as window ends!', mfilename);
else
	nwin = length(winstart);
end

if isempty(timestamps)
	error('%s: no timestamps!!!!!', mfilename)
end

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