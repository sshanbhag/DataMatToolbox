%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

% make local copy of Marker structure
Marker = D.Marker;

% get # of markers by looking at length of the Timestamp vector 
% could use any of the vectors
Nmarkers = length(Marker.Timestamp);

%-----------------------------------------------------------------------------
% build a global stimulus tag cell array, Tags
% and an vector of strings that will be compared
%-----------------------------------------------------------------------------
% preallocate Tags cell array to hold data from all Marker stimulus "tags"
% (a.k.a. stimulus parameters) and TagStrings{} cell vector to store
% strings made up from each Marker's Tags (for use in searching, see below)
Tags = cell(Nmarkers, 2*NTAGS_PER_CHANNEL);
TagStrings = cell(Nmarkers, 1);

% ctxt will be used to select Right and Left channels in the
% STIMULUS_SEARCH_TAGS{} strings
ctxt = 'RL';

% loop through Markers
for m = 1:Nmarkers
	% loop through number of tags
	for t = 1:NTAGS_PER_CHANNEL
		% loop through channels
		for c = 1:N_CHANNELS
			% offset to handle R vs L channel data
			offset = (c-1)*NTAGS_PER_CHANNEL;

			% need to check if each specific tag is in a cell or vector in order
			% to address the tag appropriately (using {} for cell, () for
			% vector)
			if iscell(Marker.([STIMULUS_SEARCH_TAGS{t} ctxt(c)]))
				Tags{m, t + offset} = Marker.([STIMULUS_SEARCH_TAGS{t} ctxt(c)]){m};
			else
				Tags{m, t + offset} = Marker.([STIMULUS_SEARCH_TAGS{t} ctxt(c)])(m);
			end
			
		end % END OF c (channel) LOOP
	end	% END OF t (tag) LOOP
end	% END OF m (marker) LOOP

%-----------------------------------------------------------------------------
% Now we have a cell array that has a mix of text and numeric fields.
% To enable a straightforward search for unique stimuli:
%	(1)	use cell2str to convert each row of Tags{} to a string that
%			stored in TagStrings{}
%	(2)	use findUniqueText function to locate unique stimuli
%-----------------------------------------------------------------------------

% perform conversion to strings.
TagStrings = cell2str(Tags);

% find unique strings
% uniqueIndices holds the indices (of Markers) for each group of common 
% stimulus parameters.
[uniqueText, uniqueIndices, Nunique] = findUniqueText(TagStrings)



%-----------------------------------------------------------------------------
% Determine Stimulus Type
%-----------------------------------------------------------------------------
% preallocate some storage structures
BG = struct('indices', [], 'count', 0);
Stimulus = repmat( cell2struct(STIMULUS_STRUCT_FIELDS, STIMULUS_STRUCT_FIELDS, 2), Nunique, 1);

for n = 1:Nunique
	
	% extract type of stimulus from Marker Data
	% we'll use the first instance of the stimulus, knowing that
	% if the sorting method to determine unique stimuli from above
	% messes up, this will also be fatally flawed...  
	
	% get index to first instance in Marker vectors
	mIndex = uniqueIndices{n}(1);
	rtype = Marker.StimulusTypeR{1, 1};
	ltype = Marker.StimulusTypeL{1, 1};
	
	% first, determine if this is a background trial (no sound played)
	% as well as the channel for the sound stimulus (left, right, both)
	if strcmp(rtype, 'NO_SOUND') && strcmp(ltype, 'NO_SOUND')
		% if both of the type tags are 'NO_SOUND', then 
		% this is a background trial
		% save Marker indices in the BG (Background) struct and increment counter
		BG.count = BG.count + 1;
		BG.indices{BG.count} = unique_indices{n};
		Stimulus(n).Type = {'BACKGROUND'};		
		Stimulus(n).Channel = 'B';
		
	elseif strcmp(rtype, 'NO_SOUND') && ~strcmp(ltype, 'NO_SOUND')
		% sound is on Left Channel
		Stimulus(n).Channel = 'L';
		Stimulus(n).Type = {ltype};
		
	elseif ~strcmp(rtype, 'NO_SOUND') && strcmp(ltype, 'NO_SOUND')
		% sound is on Right Channel
		Stimulus(n).Channel = 'R';
		Stimulus(n).Type = {rtype};
	else
		% Both channels
		Stimulus(n).Channel = 'B';
		Stimulus(n).Type = {rtype ltype};
	end
	
	Stimulus(n).Indices = uniqueIndices{n};
	Stimulus(n).tagstring = uniqueText{n};
	
	% Now copy over the Marker tags
	
end

%-----------------------------------------------------------------------------
% Now, for each stimulus, figure out what variables
%-----------------------------------------------------------------------------
