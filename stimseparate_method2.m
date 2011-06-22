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
% preallocate Tags cell array
Tags = cell(Nmarkers, 2*NTAGS_PER_CHANNEL);
TagStrings = cell(Nmarkers, 1);


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
			% to address the tag appropriately (using {} for cell, () for vector)
		
			if iscell(Marker.([STIMULUS_TAGS{t} ctxt(c)]))
				Tags{m, t + offset} = Marker.([STIMULUS_TAGS{t} ctxt(c)]){m};
			else
				Tags{m, t + offset} = Marker.([STIMULUS_TAGS{t} ctxt(c)])(m);
			end
		end
	end
end


[tmp, Ncols] = size(Tags);

% Now we have a cell array that has a mix of text and numeric fields.
% To enable a straightforward search for unique stimuli:
%	(1)	use cell2str to convert each row of Tags{} to a string that
%			stored in TagStrings{}
%	(2)	use findUniqueText function to locate unique stimuli

% perform conversion to strings.
TagStrings = cell2str(Tags);

% find unique strings
% uniqueIndices holds the indices (of Markers) for each group of common 
% stimulus parameters.
[uniqueText, uniqueIndices, Nunique] = findUniqueText(TagStrings)



%  Now  pull together  all the stimuli that differ only by attenuation value
for n = 1:Nunique
end


