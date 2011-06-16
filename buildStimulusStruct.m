%------------------------------------------------------------------------
% [UnitData, ProbeData, errFlg] = parseDataWaveProbes(ProbeData, Marker)
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%
%
%------------------------------------------------------------------------
% See:  
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 13 June, 2011 (SJS)
% 	- uses code snipped from loadDWPLXFileAsRaster.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

load test.mat -MAT

% First, need to figure out what the different stimulus conditions are
% 
% StimulusID
% 
% 	type
% 	test
% 	channel
% 	
% 
% STIMULUS_TYPES = { ...
% 	'NOISE', ...
% 	'WAVFILE', ...
% 	'TONE' ...
% 	'NO_SOUND', ...
% 	'UNKNOWN' ...
% };

% get # of markers
Nmarkers = length(D.Marker.Timestamp);

% build stimulus type cell array
stimulus_tags = [D.Marker.StimulusTypeR D.Marker.StimulusTypeL];

stimulus_tags{3, 2} = 'TONE';

% generate stimulus id # array - this is because the MATLAB unique()
% function cannot deal with finding unique rows of elements in 
% cell arrays.  so here, loop through the text tags and convert them
% to integers corresponding to their location in the STIMULUS_TYPES
% cell vector.
stimulus_ids = zeros(size(stimulus_tags));
for n = 1:prod(size(stimulus_tags))
	stimulus_ids(n) = idtag2idnum(stimulus_tags(n), STIMULUS_TYPES);
end

% now find the first unique rows of stimuli which then correspond to 
% unique L and R channel stimulus combinations
[unique_ids, locs] = unique(stimulus_ids, 'rows', 'first')
	
Nunique = length(unique_ids);

if isempty(Nunique)
	error('%s: no unique stimuli found', mfilename);
end

% Convert back to tags
unique_tags = cell(size(unique_ids));
for n = 1:prod(size(unique_ids))
	unique_tags{n} = idnum2idtag(unique_ids(n), STIMULUS_TYPES);
end

% Now find the locations of the unique stimulus combinations in the 
% stimulus_ids
unique_indices = cell(Nunique, 1);
for n = 1:Nunique
	hits = zeros(Nmarkers, 1);
	for m = 1:Nmarkers
		hits(m) = all(stimulus_ids(m, :) == unique_ids(n, :));
	end
	unique_indices{n} = find(hits);
end

%-----------------------------------------------------------------------------
% Now, we have found the unique STIMULUS TYPES in the stimuli presented
% for each channel.  However, it is still not known if there are differences 
% within the stimulus types; e.g., attenuation values may differ between 
% stimuli, frequency bandwith might vary, etc.  This next chunk of code will
% attempt to sort this out
%-----------------------------------------------------------------------------
Marker = D.Marker;

NTAGS_PER_CHANNEL = length(STIMULUS_TAGS);

% loop through unique indices
for n = 1:Nunique
	% loop through the indices for this set of stimuli
	for m = 1:length(unique_indices{n})
		u = unique_indices{n}(m);
		for t = 1:NTAGS_PER_CHANNEL
			if iscell(Marker.([STIMULUS_TAGS{t} 'R']))
				T{m, t} = Marker.([STIMULUS_TAGS{t} 'R']){u};
			else
				T{m, t} = Marker.([STIMULUS_TAGS{t} 'R'])(u);
			end
			
			if iscell(Marker.([STIMULUS_TAGS{t} 'L']))
				T{m, t+length(STIMULUS_TAGS)} = Marker.([STIMULUS_TAGS{t} 'L']){u};
			else
				T{m, t+length(STIMULUS_TAGS)} = Marker.([STIMULUS_TAGS{t} 'L'])(u);
			end
		end
	end
	StimTags{n} = T;
end


% build a global stimulus tag array
for m = 1:Nmarkers
	for t = 1:NTAGS_PER_CHANNEL
		if iscell(Marker.([STIMULUS_TAGS{t} 'R']))
			G{m, t} = Marker.([STIMULUS_TAGS{t} 'R']){m};
		else
			G{m, t} = Marker.([STIMULUS_TAGS{t} 'R'])(m);
		end

		if iscell(Marker.([STIMULUS_TAGS{t} 'L']))
			G{m, t+length(STIMULUS_TAGS)} = Marker.([STIMULUS_TAGS{t} 'L']){m};
		else
			G{m, t+length(STIMULUS_TAGS)} = Marker.([STIMULUS_TAGS{t} 'L'])(m);
		end
	end
end


% preallocate some storage structures
BG = struct('indices', [], 'count', 0);
Var = repmat(	...
					struct(	'ncols',			0,		...
								'cols',			[],	...
								'uniquevals',	cell(1, 1)	...
							), ...
					Nunique, ...
					1	...
				);
			
% initialize stimulus counter and channel cell array
StimCount = 0;
channel = cell(Nunique,	1);
stimlist = cell(1, 1);

% loop through unique stimulus types
for n = 1:Nunique
	T = StimTags{n};
	[nrows, ncols] = size(T);

	if strcmp(T{n, 1}, 'NO_SOUND') && ~strcmp(T{n, NTAGS_PER_CHANNEL+1}, 'NO_SOUND')
		% if both of the type tags are 'NO_SOUND', then 
		% this is a background trial
		% save Marker indices in the BG (Background) struct and increment counter
		BG.indices = unique_indices{n};
		BG.count = BG.count + 1;
		% set tcols to empty so that later processing stages can skip further
		% analyses of the tag column values
		tcols = [];
		channel{n} = 'B';
	% if either of the type tags are 'NO_SOUND', then that section can be ignored
	elseif strcmp(T{n, 1}, 'NO_SOUND') && ~strcmp(T{n, NTAGS_PER_CHANNEL+1}, 'NO_SOUND')
		% set tcols to Tags corresponding to L channel only
		tcols = (NTAGS_PER_CHANNEL + 2):(2*NTAGS_PER_CHANNEL);
		channel{n} = 'L';
	elseif ~strcmp(T{n, 1}, 'NO_SOUND') && strcmp(T{n, NTAGS_PER_CHANNEL+1}, 'NO_SOUND')
		% set tcols to Tags corresponding to R channel only
		tcols = (2:NTAGS_PER_CHANNEL);
		channel{n} = 'R';
	else
		% set tcols to Tags corresponding to both L and R channels
		tcols = [2:NTAGS_PER_CHANNEL (NTAGS_PER_CHANNEL + 2):(2*NTAGS_PER_CHANNEL)];
		channel{n} = 'B';
	end


	if ~isempty(tcols)
		if nrows == 1
			% this is by necessity its own stimulus type
			StimCount = StimCount + 1;
			stimlist{StimCount} = unique_indices{n};
		else
			% now look for differences in different tags
			
			% tag loop
			for col = tcols
				testvals = cell2mat(T(:, col));
				unique_testvals = unique(testvals, 'first');
				if length(unique_testvals) > 1
					Var(n).ncols = Var(n).ncols + 1;
					Var(n).cols(Var(n).ncols) = col;
					Var(n).uniquevals{Var(n).ncols} = unique_testvals;
				end
			end
		end
	end
end


%-----------------------------------------------------------------------------
% Here's where some assumptions about the stimuli need to be made.
%
%
%-----------------------------------------------------------------------------
	

% 
% 
% 
% %-----------------------------------------------------------------------------
% % build Stimulus structure
% % Each element of the Stimulus struct array contains information for 
% % a specific stimulus (either .wav file or frequency or bandwidth)
% % 
% % 	Field			Description
% % 	name				.wav file name for stimulus, or frequency for freq-resp area
% % 	indices			locations in Marker array for this stimulus
% % 	attenvals		stimulus attenuation values used for this stimulus
% % 	sweepindex		cell array of vectors that hold indices for arrays in
% % 						the Marker structure where this stimulus was played. Has
% % 						length(attenvals) elements, with each cell element corresponding
% % 						to the respective attenuation value.
% % 						Therefore, for attenvals(2), sweepindex{2} will hold the
% % 						indices for Marker  arrays corresponding to Stimulus.name
% % 						being played at attenvals(2) attenuation
% %-----------------------------------------------------------------------------
% Stimulus = struct('name', uniqueText, 'indices', uniqueIndices);
% 
% %-----------------------------------------------------------------------------
% % Find attenuation values
% %-----------------------------------------------------------------------------
% 
% % initialize attenuation arrays to store attenuation values for 
% % each of the unique stimulus strings
% attvals = cell(NuniqueStim, 1);
% attvalsUnique = cell(NuniqueStim, 1);
% 
% % Now, for each unique stimulus, find the attenuation levels
% for s = 1:NuniqueStim
% 	% get # of sweeps for this stimulus from length of the indices vector
% 	Stimulus(s).nsweeps = length(Stimulus(s).indices);
% 	% save timestamps for this stimulus
% 	Stimulus(s).tstamps = Marker.tstamps(Stimulus(s).indices);
% 	
% 	% retrieve the attenuation values used with this stimulus
% 	attvals{s} = Marker.value(Stimulus(s).indices);
% 	% determine unique attenuation values for this stimulus
% 	attvalsUnique{s} = remdupl(attvals{s}, 1e-10, 0);
% 	% store in stimulus struct array
% 	Stimulus(s).attenvals = attvalsUnique{s};
% 	Stimulus(s).n_attenvals = length(attvalsUnique{s});
% 
% 	% Find sweep index for this stimulus & for each atten value
% 	% we'll first need the locations of the stimuli
% 	list1 = strcmp(Stimulus(s).name, Marker.text);
% 	
% 	% pre-allocate some variables
% 	Stimulus(s).sweepindex = cell(Stimulus(s).n_attenvals, 1);
% 	Stimulus(s).sweeps_per_atten = zeros(Stimulus(s).n_attenvals, 1);
% 	Stimulus(s).tstamps_per_atten = cell(Stimulus(s).n_attenvals, 1);
% 	Stimulus(s).sweepstart_t = cell(Stimulus(s).n_attenvals, 1);
% 	Stimulus(s).sweepend_t = cell(Stimulus(s).n_attenvals, 1);
% 	
% 	% loop throught the # of attenuation values for this stimulus
% 	for n = 1:Stimulus(s).n_attenvals
% 		% find current attenuation values in the values list
% 		list2 = (Stimulus(s).attenvals(n) == Marker.value);
% 		% AND the two lists
% 		list1AND2 = list1 & list2;
% 		% get the indices for this list
% 		Stimulus(s).sweepindex{n} = find(list1AND2);
% 		Stimulus(s).sweeps_per_atten(n) = length(Stimulus(s).sweepindex{n});
% 		% store indices in global list
% 		attvalsIndices{s, n} = Stimulus(s).sweepindex{n};
% 		% save timestamps for this stimulus/atten combination
% 		Stimulus(s).tstamps_per_atten{n} = Marker.tstamps(Stimulus(s).sweepindex{n});
% 		
% 		% compute start and end time for this sweep - this will be used
% 		% to locate the spike times that occurred during this stimulus
% 		% presentation window
% 		%---------------------------------------------------------------------
% 		% first, get the timestamps for this stimulus/atten combination - 
% 		% this defines the start of the sweep time, store in sweepstart_t cell
% 		
% 		Stimulus(s).sweepstart_t{n} = Stimulus(s).tstamps_per_atten{n};
% 		Stimulus(s).sweepend_t{n} = 0 * Stimulus(s).sweepstart_t{n};
% 		for m = 1:Stimulus(s).sweeps_per_atten(n)
% 			% check if this is the last sweep in the experiment
% 			if Stimulus(s).sweepstart_t{n}(m) == Marker.tstamps(end)
% 				% if so, use some default number to it to calculate the end time
% 				% (in this case, use the cumulative sum of the marker times - this
% 				% will ensure that scale of the default end time is on the order of
% 				% the tstamps units
% 				Stimulus(s).sweepend_t{n}(m) = sum(Marker.tstamps);
% 			else
% 				% if not, use the time of the next timestamp in Marker.tstamps,
% 				% which is given by adding 1 to the index for this sweep
% 				Stimulus(s).sweepend_t{n}(m) = Marker.tstamps(Stimulus(s).sweepindex{n}(m) + 1);
% 			end
% 		end
% 		
% 	end
% end
% 
% % get a list of all attenuation values used
% attenList = remdupl(cell2mat(attvals), 1e-10, 0);
% 
% % get # of sweeps per stimulus
% sweepCount = cell(NuniqueStim, 1);
% for s = 1:NuniqueStim
% 	sweepCount{s} = Stimulus(s).sweeps_per_atten;
% 	nsweeps(s) = max(sweepCount{s});
% end
% maxSweeps = max(nsweeps);
% 
% 
% % now, get the spikes for each stimulus X atten combination
% for s = 1:NuniqueStim
% 	% allocate the cell array to store valid spike times for each unit
% 	Stimulus(s).spikes = cell(plxdata.total_noof_units, Stimulus(s).n_attenvals);
% 	
% 	% loop through the units
% 	for u = 1:plxdata.total_noof_units
% 
% 		% for clarity, get the spikes for this unit
% 		unit_spikes = plxdata.ts{u};
% 
% 		% loop through attenuation values
% 		for n = 1:Stimulus(s).n_attenvals
% 			
% 			% create tmp cell vector to hold spiketimes for each sweep
% 			tmpspikes = cell(Stimulus(s).sweeps_per_atten(n), 1);
% 						
% 			% loop through sweeps
% 			for p = 1:Stimulus(s).sweeps_per_atten(n)
% 				% retrieve the spikes that are valid, using the sweepstart_t and 
% 				% sweepend_t values for this stimulus/attenuation combination
% 
% 				% spiketimes greater than sweep start
% 				above_start = (unit_spikes > Stimulus(s).sweepstart_t{n}(p));
% 				% spiketime less than sweep end
% 				below_end	= (unit_spikes < Stimulus(s).sweepend_t{n}(p));
% 				% AND the two lists
% 				valid_times_vector = above_start & below_end;
% 				% get the indices and corresponding times for this list
% 				valid_index = find(valid_times_vector);
% 				valid_times = unit_spikes(valid_index);
% 				% store the values
% 				if ~isempty(valid_index)
% 					tmpspikes{p} = valid_times;
% 				else
% 					tmpspikes{p} = [];
% 				end
% 			end
% 			
% 			% store tmpspikes in the Stimulus cell array
% 			Stimulus(s).spikes{u, n} = tmpspikes;
% 		end
% 	end
% end