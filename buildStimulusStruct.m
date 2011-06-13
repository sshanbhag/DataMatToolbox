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



% First, need to figure out what the different stimulus conditions are
% 














%-----------------------------------------------------------------------------
% build Stimulus structure
% Each element of the Stimulus struct array contains information for 
% a specific stimulus (either .wav file or frequency or bandwidth)
% 
% 	Field			Description
% 	name				.wav file name for stimulus, or frequency for freq-resp area
% 	indices			locations in Marker array for this stimulus
% 	attenvals		stimulus attenuation values used for this stimulus
% 	sweepindex		cell array of vectors that hold indices for arrays in
% 						the Marker structure where this stimulus was played. Has
% 						length(attenvals) elements, with each cell element corresponding
% 						to the respective attenuation value.
% 						Therefore, for attenvals(2), sweepindex{2} will hold the
% 						indices for Marker  arrays corresponding to Stimulus.name
% 						being played at attenvals(2) attenuation
%-----------------------------------------------------------------------------
Stimulus = struct('name', uniqueText, 'indices', uniqueIndices);

%-----------------------------------------------------------------------------
% Find attenuation values
%-----------------------------------------------------------------------------

% initialize attenuation arrays to store attenuation values for 
% each of the unique stimulus strings
attvals = cell(NuniqueStim, 1);
attvalsUnique = cell(NuniqueStim, 1);

% Now, for each unique stimulus, find the attenuation levels
for s = 1:NuniqueStim
	% get # of sweeps for this stimulus from length of the indices vector
	Stimulus(s).nsweeps = length(Stimulus(s).indices);
	% save timestamps for this stimulus
	Stimulus(s).tstamps = Marker.tstamps(Stimulus(s).indices);
	
	% retrieve the attenuation values used with this stimulus
	attvals{s} = Marker.value(Stimulus(s).indices);
	% determine unique attenuation values for this stimulus
	attvalsUnique{s} = remdupl(attvals{s}, 1e-10, 0);
	% store in stimulus struct array
	Stimulus(s).attenvals = attvalsUnique{s};
	Stimulus(s).n_attenvals = length(attvalsUnique{s});

	% Find sweep index for this stimulus & for each atten value
	% we'll first need the locations of the stimuli
	list1 = strcmp(Stimulus(s).name, Marker.text);
	
	% pre-allocate some variables
	Stimulus(s).sweepindex = cell(Stimulus(s).n_attenvals, 1);
	Stimulus(s).sweeps_per_atten = zeros(Stimulus(s).n_attenvals, 1);
	Stimulus(s).tstamps_per_atten = cell(Stimulus(s).n_attenvals, 1);
	Stimulus(s).sweepstart_t = cell(Stimulus(s).n_attenvals, 1);
	Stimulus(s).sweepend_t = cell(Stimulus(s).n_attenvals, 1);
	
	% loop throught the # of attenuation values for this stimulus
	for n = 1:Stimulus(s).n_attenvals
		% find current attenuation values in the values list
		list2 = (Stimulus(s).attenvals(n) == Marker.value);
		% AND the two lists
		list1AND2 = list1 & list2;
		% get the indices for this list
		Stimulus(s).sweepindex{n} = find(list1AND2);
		Stimulus(s).sweeps_per_atten(n) = length(Stimulus(s).sweepindex{n});
		% store indices in global list
		attvalsIndices{s, n} = Stimulus(s).sweepindex{n};
		% save timestamps for this stimulus/atten combination
		Stimulus(s).tstamps_per_atten{n} = Marker.tstamps(Stimulus(s).sweepindex{n});
		
		% compute start and end time for this sweep - this will be used
		% to locate the spike times that occurred during this stimulus
		% presentation window
		%---------------------------------------------------------------------
		% first, get the timestamps for this stimulus/atten combination - 
		% this defines the start of the sweep time, store in sweepstart_t cell
		
		Stimulus(s).sweepstart_t{n} = Stimulus(s).tstamps_per_atten{n};
		Stimulus(s).sweepend_t{n} = 0 * Stimulus(s).sweepstart_t{n};
		for m = 1:Stimulus(s).sweeps_per_atten(n)
			% check if this is the last sweep in the experiment
			if Stimulus(s).sweepstart_t{n}(m) == Marker.tstamps(end)
				% if so, use some default number to it to calculate the end time
				% (in this case, use the cumulative sum of the marker times - this
				% will ensure that scale of the default end time is on the order of
				% the tstamps units
				Stimulus(s).sweepend_t{n}(m) = sum(Marker.tstamps);
			else
				% if not, use the time of the next timestamp in Marker.tstamps,
				% which is given by adding 1 to the index for this sweep
				Stimulus(s).sweepend_t{n}(m) = Marker.tstamps(Stimulus(s).sweepindex{n}(m) + 1);
			end
		end
		
	end
end

% get a list of all attenuation values used
attenList = remdupl(cell2mat(attvals), 1e-10, 0);

% get # of sweeps per stimulus
sweepCount = cell(NuniqueStim, 1);
for s = 1:NuniqueStim
	sweepCount{s} = Stimulus(s).sweeps_per_atten;
	nsweeps(s) = max(sweepCount{s});
end
maxSweeps = max(nsweeps);


% now, get the spikes for each stimulus X atten combination
for s = 1:NuniqueStim
	% allocate the cell array to store valid spike times for each unit
	Stimulus(s).spikes = cell(plxdata.total_noof_units, Stimulus(s).n_attenvals);
	
	% loop through the units
	for u = 1:plxdata.total_noof_units

		% for clarity, get the spikes for this unit
		unit_spikes = plxdata.ts{u};

		% loop through attenuation values
		for n = 1:Stimulus(s).n_attenvals
			
			% create tmp cell vector to hold spiketimes for each sweep
			tmpspikes = cell(Stimulus(s).sweeps_per_atten(n), 1);
						
			% loop through sweeps
			for p = 1:Stimulus(s).sweeps_per_atten(n)
				% retrieve the spikes that are valid, using the sweepstart_t and 
				% sweepend_t values for this stimulus/attenuation combination

				% spiketimes greater than sweep start
				above_start = (unit_spikes > Stimulus(s).sweepstart_t{n}(p));
				% spiketime less than sweep end
				below_end	= (unit_spikes < Stimulus(s).sweepend_t{n}(p));
				% AND the two lists
				valid_times_vector = above_start & below_end;
				% get the indices and corresponding times for this list
				valid_index = find(valid_times_vector);
				valid_times = unit_spikes(valid_index);
				% store the values
				if ~isempty(valid_index)
					tmpspikes{p} = valid_times;
				else
					tmpspikes{p} = [];
				end
			end
			
			% store tmpspikes in the Stimulus cell array
			Stimulus(s).spikes{u, n} = tmpspikes;
		end
	end
end