
for n = 1:d.Stimuli.N
	fn{n} = d.Stimuli.S{n, 2}.Filename;
end

runFlag = 1;
nUnique = 0;

% search all indices initially
searchIndices = 1:d.Stimuli.N;

while runFlag
	stim = searchIndices(1);
	fprintf('%d:\t%s\t%s\n', stim, d.Stimuli.S{stim, 2}.Filename, d.Stimuli.Channel(stim))

	% check channel for current stimulus
	c = d.Stimuli.Channel(stim);
	
	% compare current stim to other stims
	if strcmpi(c, 'L') || strcmpi(c, 'B')
		% if channel is L or Both, compare left, store matches in lcomp
		[lcomp, llist] = d.Stimuli.S{stim, 1}.match(d.Stimuli.S(searchIndices, 1));
	end
	if strcmpi(c, 'R') || strcmpi(c, 'B')
		% if channel is R or Both, compare right, store matches in rcomp
		[rcomp, rlist] = d.Stimuli.S{stim, 2}.match(d.Stimuli.S(searchIndices, 2));
	end
	if strcmpi(c, 'B')
		% if channel is Both, AND the lcomp and rcomp
		comp = lcomp & rcomp;
	else
		% otherwise, use appropriate channels results
		if strcmpi(c, 'L')
			comp = lcomp;
		else
			comp = rcomp;
		end
	end
	
	% store unique indices
	nUnique = nUnique + 1;
	group{nUnique} = find(comp); %#ok<SAGROW>
	% eliminate them from the list to search
	searchIndices = searchIndices(~logical(comp));

	if (stim == d.Stimuli.N) || isempty(searchIndices)
		% stop
		runFlag = 0;
	end
	
end

