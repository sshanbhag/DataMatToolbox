%------------------------------------------------------------------------
% plotRasters
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
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%-----------------------------------------------------------
% load data
%-----------------------------------------------------------
load loadDWfile_debug.mat


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
	% loop through frequencies (for tone stimuli)
	for stimindex = 1:Nstimuli
		% loop through different attenuation values
		for attenindex = 1:length(StimList(stimindex).RAttenVals)
			% clear the tmpcell variable to avoid confusing results
			clear tmpcell;
			
			% loop through individual sweeps for this attenuation setting
			for sweepindex = 1:length(StimList(stimindex).RAttenIndices{attenindex})
				spikecol = StimList(stimindex).RAttenIndices{attenindex}(sweepindex);
				tmpcell{sweepindex} = StimList(stimindex).Spiketimes{unitindex, spikecol};
				% subtract off Spikestart(sweepindex) to get times relative to 
				% start of sweep and then convert from usec to msec
				if ~isempty(tmpcell{sweepindex})
					tmpcell{sweepindex} = tmpcell{sweepindex} -  StimList(stimindex).Sweepstart(spikecol);
					tmpcell{sweepindex} = 0.001 * tmpcell{sweepindex};
				end
			end
			Spikes{unitindex, stimindex, attenindex} = tmpcell;
		end	% end of ATTENINDEX
	end	% end of FREQINDEX
end	% end of UNITINDEX


% plot a raster for each unit
H = figure(1);

for unit = 1:D.Info.Nunits

	for var = 1:Nstimuli
	
		for atten = 1:length(StimList(stimindex).RAttenVals)
			rasterplot(Spikes{unit, stimindex, atten}, [0 1000], H);
	
			titlestr = sprintf('Unit %d: Atten = %d Stim = %.2f', ...
							unit, StimList(stimindex).RAttenVals(atten), StimList(stimindex).Var(1).values(1));
			title(titlestr)
	
			pause
		end
	end
end

		

