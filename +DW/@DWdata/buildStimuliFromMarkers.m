function Stimuli = buildStimuliFromMarkers(obj, varargin)
%------------------------------------------------------------------------
% Stimuli = buildStimuliFromMarkers(obj, varargin)
%------------------------------------------------------------------------
% DataMat Toolbox
% DW Package
% Method Definition
%------------------------------------------------------------------------
%
% Takes DataMat marker information structure and splits up data by 
% stimulus type, storing relevant Marker list in the MarkerList object
%
%------------------------------------------------------------------------
% Input Arguments:
%
% Output Arguments:
%
% 	Stimuli			Structure array of data, sorted by stimulus type
% 	 Fields:
% 		Type							'TONE', 'NOISE', 'WAVFILE'
% 		Channel						'R', 'L', 'B'
% 		Indices						array of indices into Marker arrays (above)
% 		Nreps							# of times this stimulus was presented 
% 		Var							varying variables for this stimulus
% 		Nsweeps						# of times this stimulus was presented
% 		Sweepstart					list of start times for this stimulus
% 		Sweepend						list of end times for this stimulus
% 		LAttenVals					Left channel attenuation values
% 		LAttenIndices				Indices for respective attenuation sweeps
% 		RAttenVals					Right channel attenuation values
% 		RAttenIndices				indices for R attenuaton sweeps
% 		Spiketimes					{# units, Nsweeps} cell array of unit spike
%										timestamps
%
% 		Timestamp					first occurance of this Stimulus
% 		Tagstring
%
%------------------------------------------------------------------------
% See also: loadDWfile, buildSpikes 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 4 June, 2011 (SJS)
% 	- modified from buildStimulusStruct.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

% make local copy of Markers object array
Marker = obj.Markers;

%-----------------------------------------------------------------------------
% build a global stimulus tag cell array, Tags
% and an vector of strings that will be compared
%-----------------------------------------------------------------------------
% preallocate Tags cell array to hold data from all Marker stimulus "tags"
% (a.k.a. stimulus parameters) and TagStrings{} cell vector to store
% strings made up from each Marker's Tags (for use in searching, see below)
Tags = cell(obj.Nmarkers, 2*NTAGS_PER_CHANNEL);
TagStrings = cell(obj.Nmarkers, 1);

% ctxt will be used to select Right and Left channels in the
% STIMULUS_SEARCH_TAGS{} strings
ctxt = 'RL';

% loop through Markers
for m = 1:obj.Nmarkers
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

return

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
[uniqueText, uniqueIndices, Nunique] = findUniqueText(TagStrings);


%-----------------------------------------------------------------------------
% Determine Stimulus Type
%-----------------------------------------------------------------------------
% preallocate some storage structures
BG = struct('indices', [], 'count', 0);
Stimulus = repmat( cell2struct(STIMULUS_STRUCT_FIELDS, STIMULUS_STRUCT_FIELDS, 2), Nunique, 1);
Nstimuli = length(Stimulus);

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
	Stimulus(n).Tagstring = uniqueText{n};
	Stimulus(n).Nreps = length(uniqueIndices{n});
	
	% Now copy over the Marker tags and values
	for f = 1:length(MARKER_TAGS)
		if iscell(Marker.(MARKER_TAGS{f}))
			Stimulus(n).(MARKER_TAGS{f}) = ...
										Marker.(MARKER_TAGS{f}){uniqueIndices{n}};
		else
			Stimulus(n).(MARKER_TAGS{f}) = ...
										Marker.(MARKER_TAGS{f})(uniqueIndices{n});
		end
	end
	
end

%-----------------------------------------------------------------------------
% this is a bit of a kludge - the Timestamp field in Stimulus only has the 
% first timestamp, whereas we need all timestamps for each stimulus type in
% 
% This information will be stored in the Sweepstart(1:Nreps) and
% Sweepend(1:Nreps) arrays.  Sweepstart corresponds to the sweep timestamp,
% and Sweepend is the timestamp for the next sweep.
%
% for the final sweep in the file, Sweepend isn't truly accurate, but is 
% set to that sweep's start timestamp + max interval between all previous sweeps. 
%-----------------------------------------------------------------------------

% Marker.Timestamp is in string format, so convert to double
marker_timestamp_vals = str2double(Marker.Timestamp);
% find max difference between timestamps
max_timestamp_difference = max(diff(marker_timestamp_vals));
if max_timestamp_difference <= 0
	error('%s: strange value for max_timestamp_difference %f', ...
						mfilename, max_timestamp_difference);
end

% then get the timestamps for this Stimulus from the Marker structure
for s = 1:Nstimuli
	Stimulus(s).Sweepstart = marker_timestamp_vals(Stimulus(s).Indices);
	Stimulus(s).Sweepend = zeros(size(Stimulus(s).Sweepstart));
	% Perform a switcheroo on the Timestamp field.
	Stimulus(s).FirstTimestamp = str2double(Stimulus(s).Timestamp);
	Stimulus(s).Timestamp = Stimulus(s).Sweepstart;
	
	% and determine end time for sweep
	for n = 1:length(Stimulus(s).Indices)
		start_index = Stimulus(s).Indices(n);
		if start_index ~= Marker.Nmarkers
			% Set the Sweep end time to the Timestamp of the next marker
			Stimulus(s).Sweepend(n) = marker_timestamp_vals(start_index + 1);
			
		else
			% if the index for this stimulus (into the Markers arrays) is
			% equal to the total # of markers, then this is the last marker and
			% the Sweepend time will need to be computed
			Stimulus(s).Sweepend(n) = marker_timestamp_vals(start_index) + max_timestamp_difference;

		end
		
	end	% end of n loop
end	% end of S loop

%-----------------------------------------------------------------------------
% another kludge - find the time for prior and post stimulus for each stimulus
%-----------------------------------------------------------------------------

% loop through stimuli
for s = 1:Nstimuli
	% compute pre indices by subtracting 1 from indices for this stimulus
	preIndices = uniqueIndices{s} - 1;
	
	% valid preIndices values are those that are greater than 1
	% invalid values are less than 1
	valid_pre = find(preIndices >= 1);
	invalid_pre = find(preIndices < 1);
	
	% preallocate preTimes
	preTimes = zeros(size(preIndices));
	% check if any invalid preTimes were found
	if any(invalid_pre)
		% if more than 1, throw an error
		if length(invalid_pre) > 1
			error('%s: strange preIndices{%d}', mfilename, invalid_pre)
		else
			% set the invalid time to max_timestamp_difference usec before first
			% stimulus marker
			min_time = min(marker_timestamp_vals) - max_timestamp_difference;
			% if greater than or equal to zero, use it, otherwise use 0
			if min_time >= 0
				preTimes(invalid_pre) = min_time;
			else
				preTimes(invalid_pre) = 0;
			end
		end
	end
	% store preTimes
	preTimes(valid_pre) = marker_timestamp_vals(preIndices(valid_pre));
	
	% add one to each of the indices for this stimulus to get the index
	% of the following stimulus onset
	postIndices = uniqueIndices{s} + 1;
	% valid postIndices values are those that are less than the total number
	% of markers; invalid values are greater than # of markers
	valid_post = find(postIndices <= length(marker_timestamp_vals));
	invalid_post = find(postIndices > length(marker_timestamp_vals));

	% pre-allocate postTimes array
	postTimes = zeros(size(postIndices));
	% check if any invalid PostTimes were found
	if any(invalid_post)
		% if more than 1, throw an error
		if length(invalid_post) > 1
			error('%s: strange postIndices{%d}', mfilename, invalid_post)
		else
			% set the invalid time to max_timestamp_difference usec after last
			% stimulus marker
			postTimes(invalid_post) = max(marker_timestamp_vals) + max_timestamp_difference;
		end
	end
	% store valid postTimes
	postTimes(valid_post) = marker_timestamp_vals(postIndices(valid_post));
	
	% assign to Stimulus struct
	Stimulus(s).PreSweep = preTimes;
	Stimulus(s).PostSweep = postTimes;
end

%-----------------------------------------------------------------------------
% Now, for each stimulus, figure out what variables are varying.  or were
% varied.  or the variety of varying variables.
%-----------------------------------------------------------------------------

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

for s = 1:Nstimuli
	% check the stimulus Type - it is stored as a string within a cell, so 
	% it must be addressed using Type{} format.  
	switch Stimulus(s).Type{1}
		case 'TONE'
			Stimulus(s).Var = assignVarTags(Stimulus(s), TONE_VAR_TAGS);
		case	'NOISE'
			Stimulus(s).Var = assignVarTags(Stimulus(s), NOISE_VAR_TAGS);
		case 	'WAVFILE'
			Stimulus(s).Var = assignVarTags(Stimulus(s), WAV_VAR_TAGS);
		otherwise
			error('%s: UNKNOWN STIMULUS TYPE: %s ', mfilename, Stimulus(s).Type{1}); 
	end
end


%------------------------------------------------------------------------
% Determine attenuation values for each stimulus
%------------------------------------------------------------------------
for s = 1:Nstimuli
	% get attenuation setting for Left channel if stimulus is played over
	% Both or Left only channel
	if strcmp(Stimulus(s).Channel, 'B') || strcmp(Stimulus(s).Channel, 'L')
		[Stimulus(s).LAttenVals, Stimulus(s).LAttenIndices, tmpn] = ...
													findUniqueValues(Stimulus(s).AttenuationL);
	else
		% otherwise, set to empty values
		Stimulus(s).LAttenVals = [];
		Stimulus(s).LAttenIndices = {};
	end
	
	% get attenuation setting for Right channel if stimulus is played over
	% Both or Right only channel
	if strcmp(Stimulus(s).Channel, 'B') || strcmp(Stimulus(s).Channel, 'R')
		[Stimulus(s).RAttenVals, Stimulus(s).RAttenIndices, tmpn] = ...
													findUniqueValues(Stimulus(s).AttenuationR);
	else
		% otherwise, set to empty values
		Stimulus(s).LAttenVals = [];
		Stimulus(s).LAttenIndices = {};
	end
end

%------------------------------------------------------------------------
% Retrieve the spikes for each stimulus presentation
%------------------------------------------------------------------------
% make local copy of UnitData
UnitData = Data.UnitData;
Nunits = Data.Info.Nunits;

% loop through stimuli
for s = 1:Nstimuli
	% allocate the cell array to store valid spike times for each unit
	Stimulus(s).Spiketimes = cell(Nunits, Stimulus(s).Nreps);

	% loop through sweeps
	for r = 1:Stimulus(s).Nreps
	
		% loop through the units
		for u = 1:Nunits
			
			% retrieve the spikes that are valid for each time period, using the 
			% sweepstart_t and sweepend_t values for each 
			% stimulus/attenuation combination

			%----------------------------------------------------------------
			% Spikes during this stimulus
			%----------------------------------------------------------------			
			% spiketimes greater than sweep start
			startT = (UnitData(u).timestamp >= Stimulus(s).Sweepstart(r));
			% spiketimes less than sweep end
			endT	= (UnitData(u).timestamp < Stimulus(s).Sweepend(r));
			% AND the two lists
			valid_times_list = startT & endT;
			% store the values if spikes were found, 
			if any(valid_times_list)
				Stimulus(s).Spiketimes{u, r} = UnitData(u).timestamp(valid_times_list);
			else
				%otherwise, set to empty array
				Stimulus(s).Spiketimes{u, r} = [];
			end
			
			%----------------------------------------------------------------
			% Spikes before this stimulus
			%----------------------------------------------------------------			
			% spiketimes greater than pre-sweep start
			startT = (UnitData(u).timestamp > Stimulus(s).PreSweep(r));
			% spiketimes less than sweep end
			endT	= (UnitData(u).timestamp < Stimulus(s).Sweepstart(r));
			% AND the two lists
			valid_times_list = startT & endT;
			% store the values if spikes were found, 
			if any(valid_times_list)
				Stimulus(s).Spiketimes_PreStim{u, r} = UnitData(u).timestamp(valid_times_list);
			else
				%otherwise, set to empty array
				Stimulus(s).Spiketimes_PreStim{u, r} = [];
			end

			%----------------------------------------------------------------
			% Spikes after this stimulus
			%----------------------------------------------------------------			
			% spiketimes greater than sweep end
			startT = (UnitData(u).timestamp > Stimulus(s).Sweepend(r));
			% spiketimes less than post-sweep end
			endT	= (UnitData(u).timestamp < Stimulus(s).PostSweep(r));
			% AND the two lists
			valid_times_list = startT & endT;
			% store the values if spikes were found, 
			if any(valid_times_list)
				Stimulus(s).Spiketimes_PostStim{u, r} = UnitData(u).timestamp(valid_times_list);
			else
				%otherwise, set to empty array
				Stimulus(s).Spiketimes_PostStim{u, r} = [];
			end			
			
		end	% end of U loop
	end	% end of R loop
end	% end of S loop

