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
	
	% Now copy over the Marker tags
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

marker_timestamp_vals = str2double(Marker.Timestamp);
max_timestamp_difference = max(diff(marker_timestamp_vals));

% then get the timestamps for this Stimulus from the Marker structure
for s = 1:Nstimuli
	Stimulus(s).Sweepstart = marker_timestamp_vals(Stimulus(s).Indices);
	Stimulus(s).Sweepend = zeros(size(Stimulus(s).Sweepstart));
	
	% and determine end time for sweep
	for n = 1:length(Stimulus(s).Indices)
		start_index = Stimulus(s).Indices(n);
		if start_index ~= Nmarkers
			% Set the Sweep end time to the Timestamp of the next marker
			Stimulus(s).Sweepend(n) = marker_timestamp_vals(start_index + 1);
			
		else
			% if the index for this stimulus (into the Markers arrays) is
			% equal to the total # of markers, then this is the last marker and
			% the Sweepend time will need to be computed
			Stimulus(s).Sweepend(n) = marker_timestamp_vals(start_index) + max_timestamp_difference;

		end
	end
			
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


% Var = repmat(	...
% 					struct(	'ncols',			0,		...
% 								'cols',			[],	...
% 								'uniquevals',	cell(1, 1),	...
% 								'channel',		[]		...
% 							), ...
% 					Nunique, ...
% 					1	...
% 				);

for s = 1:Nstimuli
	
	varindex = 0;
	
	
	Stimulus(s) = setfield(Stimulus(s), 'Var', struct('name', [], 'values', [], 'indices', []) );
	
	switch Stimulus(s).Type{1}

		case 'TONE'
			for n = 1:length(TONE_VAR_TAGS)
				
				clear tmpvar;
				
				if Stimulus(s).Channel == 'B'
					tmp = Stimulus(s).([TONE_VAR_TAGS{n} 'R']);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar(1).name = [TONE_VAR_TAGS{n} 'R'];
						tmpvar(1).values = uniqtmp;

						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar(1).indices = Stimulus(s).Indices(tmpind);
						end
					end
					
					
					tmp = Stimulus(s).([TONE_VAR_TAGS{n} 'L']);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar(2).name = [TONE_VAR_TAGS{n} 'L'];
						tmpvar(2).values = uniqtmp;

						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar(2).indices = Stimulus(s).Indices(tmpind);
						end
					end
									
				else
					tmp = Stimulus(s).([TONE_VAR_TAGS{n} Stimulus(s).Channel]);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar.name = [TONE_VAR_TAGS{n} Stimulus(s).Channel];
						tmpvar.values = uniqtmp;
						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar.indices = Stimulus(s).Indices(tmpind);
						end
					end
				end
				if exist('tmpvar')
					if ~isempty(tmpvar)
						varindex = varindex + 1;
						Stimulus(s).Var(varindex) = tmpvar;
					end
				end

			end
			
		case	'NOISE'

			for n = 1:length(NOISE_VAR_TAGS)
				
				clear tmpvar;
				
				if Stimulus(s).Channel == 'B'
					tmp = Stimulus(s).([NOISE_VAR_TAGS{n} 'R']);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar(1).name = [NOISE_VAR_TAGS{n} 'R'];
						tmpvar(1).values = uniqtmp;

						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar(1).indices = Stimulus(s).Indices(tmpind);
						end
					end
					
					
					tmp = Stimulus(s).([NOISE_VAR_TAGS{n} 'L']);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar(2).name = [NOISE_VAR_TAGS{n} 'L'];
						tmpvar(2).values = uniqtmp;

						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar(2).indices = Stimulus(s).Indices(tmpind);
						end
					end
									
				else
					tmp = Stimulus(s).([NOISE_VAR_TAGS{n} Stimulus(s).Channel]);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar.name = [NOISE_VAR_TAGS{n} Stimulus(s).Channel];
						tmpvar.values = uniqtmp;
						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar.indices = Stimulus(s).Indices(tmpind);
						end
					end
				end
				if exist('tmpvar')
					if ~isempty(tmpvar)
						varindex = varindex + 1;
						Stimulus(s).Var(varindex) = tmpvar;
					end
				end

			end
			
		case 	'WAVFILE'
			for n = 1:length(WAV_VAR_TAGS)
				
				clear tmpvar;
				
				if Stimulus(s).Channel == 'B'
					tmp = Stimulus(s).([WAV_VAR_TAGS{n} 'R']);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar(1).name = [WAV_VAR_TAGS{n} 'R'];
						tmpvar(1).values = uniqtmp;

						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar(1).indices = Stimulus(s).Indices(tmpind);
						end
					end
					
					
					tmp = Stimulus(s).([WAV_VAR_TAGS{n} 'L']);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar(2).name = [WAV_VAR_TAGS{n} 'L'];
						tmpvar(2).values = uniqtmp;

						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar(2).indices = Stimulus(s).Indices(tmpind);
						end
					end
									
				else
					tmp = Stimulus(s).([WAV_VAR_TAGS{n} Stimulus(s).Channel]);
					uniqtmp = unique(tmp);
					if ~isempty(uniqtmp)
						tmpvar.name = [WAV_VAR_TAGS{n} Stimulus(s).Channel];
						tmpvar.values = uniqtmp;
						for u = 1:length(uniqtmp)
							tmpind = find(tmp == uniqtmp(u));
							tmpvar.indices = Stimulus(s).Indices(tmpind);
						end
					end
				end
				if exist('tmpvar')
					if ~isempty(tmpvar)
						varindex = varindex + 1;
						Stimulus(s).Var(varindex) = tmpvar;
					end
				end

			end
			
			
		otherwise
			error('%s: UNKNOWN STIMULUS TYPE: %s ', mfilename, Stimulus(s).Type{1}); 
	end
end



%------------------------------------------------------------------------
% Determine attenuation values for each stimulus
%------------------------------------------------------------------------


for s = 1:Nstimuli
	if strcmp(Stimulus(s).Channel, 'B') || strcmp(Stimulus(s).Channel, 'L')
		[Stimulus(s).LAttenVals, Stimulus(s).LAttenIndices, tmpn] = ...
													findUniqueValues(Stimulus(s).AttenuationL);
	else
		Stimulus(s).LAttenVals = [];
		Stimulus(s).LAttenIndices = {};
	end
	
	if strcmp(Stimulus(s).Channel, 'B') || strcmp(Stimulus(s).Channel, 'R')
		[Stimulus(s).RAttenVals, Stimulus(s).RAttenIndices, tmpn] = ...
													findUniqueValues(Stimulus(s).AttenuationR);
	else
		Stimulus(s).LAttenVals = [];
		Stimulus(s).LAttenIndices = {};
	end
	
end


%------------------------------------------------------------------------
% Retrieve the spikes for each stimulus X atten combination
%------------------------------------------------------------------------

% make local copy of UnitData
UnitData = D.UnitData;
Nunits = length(UnitData);

for s = 1:Nstimuli
	% allocate the cell array to store valid spike times for each unit
	Stimulus(s).Spiketimes = cell(Nunits, Stimulus(s).Nreps);

	% loop through reps
	for r = 1:Stimulus(s).Nreps
	
		% loop through the units
		for u = 1:Nunits
			% retrieve the spikes that are valid, using the sweepstart_t and 
			% sweepend_t values for this stimulus/attenuation combination

			% spiketimes greater than sweep start
			above_start = (UnitData(u).timestamp >= Stimulus(s).Sweepstart(r));
			% spiketime less than sweep end
			below_end	= (UnitData(u).timestamp < Stimulus(s).Sweepend(r));
			
			% AND the two lists
			valid_times_list = above_start & below_end;
			% get the indices and corresponding times for this list
			valid_index = find(valid_times_list);

			% store the values
			if ~isempty(valid_index)
				Stimulus(s).Spiketimes{u, r} = UnitData(u).timestamp(valid_index);
			else
				Stimulus(s).Spiketimes{u, r} = [];
			end
		end
	end
end






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
