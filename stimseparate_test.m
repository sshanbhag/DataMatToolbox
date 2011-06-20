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

D.Marker.StimulusTypeR{3} = 'TONE';

% build stimulus type cell array
stimulus_tags = [D.Marker.StimulusTypeR D.Marker.StimulusTypeL];

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
	unique_indices{n}
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
	clear T
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
			G{m, t+NTAGS_PER_CHANNEL} = Marker.([STIMULUS_TAGS{t} 'L']){m};
		else
			G{m, t+NTAGS_PER_CHANNEL} = Marker.([STIMULUS_TAGS{t} 'L'])(m);
		end
	end
end


% preallocate some storage structures
BG = struct('indices', [], 'count', 0);
			
% initialize stimulus counter and channel cell array
channel = cell(Nunique,	1);

% loop through unique stimulus types
for n = 1:Nunique
	% make a local copy of the Tags matrix for this stimulus category
	T = StimTags{n};
	% and get size of it
	[nrows, ncols] = size(T);
	% extract type of stimulus for each channel
	rtype = T{1, 1};
	ltype = T{1, NTAGS_PER_CHANNEL+1};
	% store these stimulus types
	StimTypes{n} = {rtype, ltype};
	
	% first, determine if this is a background trial (no sound played)
	% as well as the channel for the sound stimulus (left, right, both)
	if strcmp(rtype, 'NO_SOUND') && strcmp(ltype, 'NO_SOUND')
		% if both of the type tags are 'NO_SOUND', then 
		% this is a background trial
		% save Marker indices in the BG (Background) struct and increment counter
		BG.count = BG.count + 1;
		BG.indices{BG.count} = unique_indices{n};
		channel{n} = 'B';
		bgFlag = 1;
	elseif strcmp(rtype, 'NO_SOUND') && ~strcmp(ltype, 'NO_SOUND')
		% sound is on Left channel
		channel{n} = 'L';
		bgFlag = 0;
	elseif ~strcmp(rtype, 'NO_SOUND') && strcmp(ltype, 'NO_SOUND')
		% sound is form Right channel
		channel{n} = 'R';
		bgFlag = 0;
	else
		% Both channels
		channel{n} = 'B';	
		bgFlag = 0;
	end		
		
		
	% if this isn't a background trial, get columns to sort for each
	% stimulus type
	if ~bgFlag
		switch channel{n}
			case 'R'
				% set tcols to Tags corresponding to R channel only
				col_offset = 0
% 				tcols{n} = (2:NTAGS_PER_CHANNEL);
			case 'L'
 				% set tcols to Tags corresponding to L channel only
				col_offset = NTAGS_PER_CHANNEL
% 				tcols{n} = (NTAGS_PER_CHANNEL + 2):(2*NTAGS_PER_CHANNEL);
			case 'B'
				% set tcols to Tags corresponding to both L and R channels
				col_offset = [0 NTAGS_PER_CHANNEL]
% 				tcols{n} = [2:NTAGS_PER_CHANNEL (NTAGS_PER_CHANNEL + 2):(2*NTAGS_PER_CHANNEL)];
		end

		tagindex = 1;
		% loop through the values in col_offset.  this will allow relevant
		% columns for particular stimuli to be retreived
		for c = col_offset
			
			% get type of stimulus for current channel
			tmptype = T{1, c + 1}
			
			switch tmptype
				case 'TONE'
					for tag = 1:length(TONE_VAR_TAGS)
						tcols{n}(tagindex) = find(strcmp(STIMULUS_TAGS, TONE_VAR_TAGS{tag})) + c;
						tagindex = tagindex + 1;
					end

				case 'NOISE'
					for tag = 1:length(NOISE_VAR_TAGS)
						tcols{n}(tagindex) = find(strcmp(STIMULUS_TAGS, NOISE_VAR_TAGS{tag})) + c;
						tagindex = tagindex + 1;
					end
					
				case 'WAV'
					for tag = 1:length(WAV_VAR_TAGS)
						tcols{n}(tagindex) = find(strcmp(STIMULUS_TAGS, WAV_VAR_TAGS{tag})) + c;
						tagindex = tagindex + 1;
					end
					
				otherwise
					tcols{n} = [];
			end
				
		end
	else
		% set tcols to empty so that later processing stages can skip further
		% analyses of the tag column values
		tcols{n} = [];	
	end
end

%-----------------------------------------------------------------------------
%
% Search for appropriately varying variables for each stimulus type
%	One might ask why this is necessary, since stimulus types have
%	already been detected.  The issue is that for some stimuli, there
%	might be characteristics that vary so that there are different "subclasses"
%	of stimuli.  E.g., the stimulus type might be WAVFILE but there might be
%	5 different wav files randomly presented - this needs to be detected
%	and accounted for.
%-----------------------------------------------------------------------------	
StimCount = 0;
stimlist = cell(1, 1);
Var = repmat(	...
					struct(	'ncols',			0,		...
								'cols',			[],	...
								'uniquevals',	cell(1, 1),	...
								'channel',		[]		...
							), ...
					Nunique, ...
					1	...
				);


% Var(1:Nunique)		struct containing information about variables
% Var(n).ncols = Var(n).ncols + 1;
% Var(n).cols(Var(n).ncols) = col;
% Var(n).uniquevals{Var(n).ncols} = unique_testvals;
% Var(n).channel = channel{n};

for n = 1:Nunique
	T = StimTags{n};
	[nrows, ncols] = size(T);

	if ~isempty(tcols{n})
		if nrows == 1
			% this is by necessity its own stimulus type
			StimCount = StimCount + 1;
			stimlist{StimCount} = unique_indices{n};
		else
			% now look for differences in different tags
			
			% tag loop
			for col = tcols{n}
				testvals = cell2mat(T(:, col));
				unique_testvals = unique(testvals, 'first');
				if length(unique_testvals) > 1
					Var(n).ncols = Var(n).ncols + 1;
					Var(n).cols(Var(n).ncols) = col;
					Var(n).uniquevals{Var(n).ncols} = unique_testvals;
					Var(n).channel = channel{n};
					Var(n).name = STIMULUS_TAGS{col};
				end
			end
			
			if Var(n).ncols == 0
				% only one configuration for this stimulus
				StimCount = StimCount + 1;
				stimlist{StimCount} = unique_indices{n};
			end
			
		end
	end
end

%-----------------------------------------------------------------------------
% collect attenuation information on stimuli in stimlist
%-----------------------------------------------------------------------------
for stimIndex = 1:StimCount
	
	testatt = Marker.stimlist{
	
	for s = 1:stimIndex
		testvals = stimlist{s


	end
end

%-----------------------------------------------------------------------------
% Here's where some assumptions about the stimuli need to be made.
%
%
%-----------------------------------------------------------------------------



for n = 1:Nunique
	if Var(n).ncols
		candidate_vars{n} = STIMULUS_TAGS{Var(n).cols};
	else
		candidate_vars{n} = [];
	end
	
	if ~isempty(candidate_vars{n})
		disp('Possible variables:')
		sprintf('%s\t', candidate_vars{n})
	end
end



