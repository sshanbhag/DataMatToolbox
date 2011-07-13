%------------------------------------------------------------------------
% plotData
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
% Created: 	4 July 2011 (SJS):
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% for test data, all stimuli are tones, so this next section
% is not strictly necessary.  however, for future use with
% intermingled stimuli, this will show how to pull out stimuli
% of specific types
%-----------------------------------------------------------
toneIndex = 0;
wavIndex = 0;
noiseIndex = 0;
unknownIndex = 0;

for s = 1:length(Stimulus)
	disp(sprintf('Stimulus(%d):\n\t type:\t %s', s, Stimulus(s).Type{1}))
	disp(sprintf('\t channel:\t %s', Stimulus(s).Channel))
	% save tone stimulus information in StimList struct array
	switch Stimulus(s).Type{1}
		case 'TONE'
			toneIndex = toneIndex + 1;
			ToneStim(toneIndex) = Stimulus(s);
			ToneFreq(toneIndex) = StimList(s).Var(1).values;
		case 'WAVFILE'
			wavIndex = wavIndex + 1;
			WavStim(wavIndex) = Stimulus(s);
		case 'NOISE'
			noiseIndex = noiseIndex + 1;
			NoiseStim(noiseIndex) = Stimulus(s);
		otherwise
			unknownIndex = unknownIndex + 1;
			UnknownStim(unknownIndex) = Stimulus(s);
	end
end

StimList = WavStim;

% get number of stimuli
Nstimuli = length(StimList);

% find ranges of attenuation
for s = 1:Nstimuli
	% R channel attenuation
	AttenR{s} = StimList(s).RAttenVals;
	AttenR_indices{s} = StimList(s).RAttenIndices;
end


%----------------------------------------------------------------------
%----------------------------------------------------------------------
% build up Spikes cell array to hold spike time information for 
%----------------------------------------------------------------------
%----------------------------------------------------------------------
% 
% Within each Stimulus (or, StimList in this case) struct element, there are 
% a few key elements for organizing the spikes.
% 
% {R/L}AttenIndices{} contains a list of indices to Stimulus vectors that correspond
% to attenuation values given in {R/L}AttenVals
% 
% for example, 
% 	
% 	The vector at Stimulus(1).RAttenIndices{1} correspondes to a list of 
% 	indices for Stimulus(1) where Stimulus(1)RAttenVals(1) was the 
% 	attenuation setting.  
% 	
% 	So, in order to get the spikestimes corresponding to Stimulus(1) and the 
% 	attenuation value of Stimulus(1).RAttenVals(1), we would use the following
% 	code:
% 	
% 	spikes = Stimulus(1).Spiketimes(<unit number>, Stimulus(1).RAttenIndices{1})
% 		(*** parentheses notation   ^     and                                   ^)
% 		
% 	Which returns a cell array of spike time stamp vectors:
% 	
% 	>> Stimulus(1).Spiketimes(1, Stimulus(1).RAttenIndices{1})
% 
% 		ans = 
% 
% 		 []    [1x15 double]    [1x14 double]
% 
% 	To obtain individual vectors of spike times, use the following notation
% 	(pulling in the second timestamp vector because the first one is empty):
% 	
% 	spikes = Stimulus(1).Spiketimes{1, Stimulus(1).RAttenIndices{1}(2)}
% 	
% 	>>	spikes = Stimulus(1).Spiketimes{1, Stimulus(1).RAttenIndices{1}(2)}
% 
% 		spikes =
% 
% 		Columns 1 through 10
% 
% 		  2868100     2869800     2870000     2870100     2870300     2870500     2870700     2870800     2871000     2871100
% 
% 		Columns 11 through 15
% 
% 		  2871200     2871300     2871500     2871600     2871700
%----------------------------------------------------------------------
% loop through units
for unitindex = 1:D.Info.Nunits
	
	% loop through stimuli (this will thus vary across frequencies or wavfiles
	for stimindex = 1:Nstimuli
		
		% loop through different attenuation values
		for attenindex = 1:length(StimList(stimindex).RAttenVals)
			
			% clear the tmpcell variable to avoid confusing results
			clear tmpcell;
			
			% loop through individual sweeps for this attenuation setting
			for sweepindex = 1:length(StimList(stimindex).RAttenIndices{attenindex})
				% To get the spikes for this sweep, we need to access the index that
				% corresponds to the proper attenuation and stimulus combination.
				% These indices are stored in 
				%	Stimulus(1... nstimuli).RattenIndices{1... # att values}(1...nsweeps)
				% 
				% retrieve this index and store it in spikecol
				spikecol = StimList(stimindex).RAttenIndices{attenindex}(sweepindex);
				
				% assign the spike times to the current tmpcell{} element
				tmpcell{sweepindex} = StimList(stimindex).Spiketimes{unitindex, spikecol};
				% to work with the spike timestamps, need to
				% (1)	subtract off Spikestart(sweepindex) to get times relative to 
				%		start of sweep
				% (2)	then convert from usec to msec
				if ~isempty(tmpcell{sweepindex})
					tmpcell{sweepindex} = 0.001*(tmpcell{sweepindex} -  StimList(stimindex).Sweepstart(spikecol));
				end
			end
			
			% assign the tmpcell spiketimes to the master Spikes cell array
			% Spikes has three dimensions: 
			%		(1) unit index
			%		(2) stimulus index (varying parameter, e.g., wav file, tone freq)
			%		(3) attenuation value
			Spikes{unitindex, stimindex, attenindex} = tmpcell;
			
		end	% end of ATTENINDEX
	end	% end of FREQINDEX
end	% end of UNITINDEX

% count number of different attenuation settings (assume same for all
% stimuli)
Natten = length(StimList(1).RAttenVals);

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% plot!
%----------------------------------------------------------------------
%----------------------------------------------------------------------

% need to get max ISI for determining plot limits

% first, get list of Marker timestamps
tmptimes = D.MarkerTimes;

% take differences
dt = diff(tmptimes);

% nonzero values
nonzero_dt = dt( dt ~= 0 );

maxtime = 0.001 * max(nonzero_dt);
mintime = 0.001 * min(nonzero_dt);




% pre-allocate plot options
plotopts = struct( ...
	'timelimits',				[0 ceil(maxtime)]			, ...
	'psth_binwidth',			5					, ...
	'raster_tickmarker',		'.'				, ...
	'raster_ticksize',		12					, ...
	'horizgap',					0.05				, ...
	'vertgap',					0.055				, ...
	'plotgap',					0.0125			, ...
	'filelabel',				D.Info.file		...
);

% create list of column labels (stimulus parameter - e.g., tone freq, wav name)
plotopts.columnlabels = cell(Nstimuli, 1);
for s = 1:Nstimuli
	if strcmp(StimList(s).Type, 'TONE')
		plotopts.columnlabels{s} = sprintf('Stim = %.2f', StimList(s).Var(1).values(1));
	elseif strcmp(StimList(s).Type, 'WAVFILE')
		tmpstr = textscan(StimList(s).Var(2).values{1}, '%s', 'Delimiter', '\\');
		plotopts.columnlabels{s} = ['Stim = ' tmpstr{1}{end}];
	else
		plotopts.columnlabels{s} = ['Stim = ' num2str(StimList(s).Var(1).values(1))];
	end
end

% create list of row labels for plots (attenuation values)
plotopts.rowlabels = cell(Natten, 1);
for n = 1:Natten
	plotopts.rowlabels{n} = sprintf('%d dB', StimList(1).RAttenVals(n));
end

% allocate temporary spike cell matrix
tmpspikes = cell(Natten, Nstimuli);

for unit = 1:D.Info.Nunits
	% for each unit, create a new page
	Hfig(unit) = figure(unit);
	
	for row = 1:Natten
		for col = 1:Nstimuli
			% assign spikes for this unit to tmpspikes - need to put different
			% attenuation levels across rows and different stimuli (tone freq,
			% wav file name, etc) across columns
			tmpspikes{row, col} = Spikes{unit, col, row};
		end
	end
	
	% write unit ID label string
	plotopts.idlabel = sprintf('Unit %d', unit);
	
	% plot a raster and psth for each unit
	[Hplots, plotopts_out] = rasterpsthmatrix(tmpspikes, plotopts);
end


