%------------------------------------------------------------------------
% 
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

Nstimuli = length(Stimulus);

%-----------------------------------------------------------
%
%-----------------------------------------------------------
for s = 1:Nstimuli
	disp(sprintf('Stimulus(%d):\n\t type:\t %s', s, Stimulus(s).Type{1}))
	disp(sprintf('\t channel:\t %s', Stimulus(s).Channel))
end


% find ranges of tone stimuli

for s = 1:Nstimuli
	if strcmp(Stimulus(s).Type, 'TONE')
		ToneFreq(s) = Stimulus(s).Var(1).values;
	end
end

% find ranges of attenuation
for s = 1:Nstimuli
	% R channel attenuation
	AttenR{s} = Stimulus(s).RAttenVals;
	AttenR_indices{s} = Stimulus(s).RAttenIndices;
end


% build up Spikes cell array to hold spike time information for 

% loop through units
for unitindex = 1:D.Info.Nunits
	% loop through frequencies (for tone stimuli)
	for freqindex = 1:Nstimuli
		% loop through different attenuation values
		for attenindex = 1:length(Stimulus(freqindex).RAttenVals)
			% loop through individual sweeps for this attenuation setting
			clear tmpcell;
			for sweepindex = 1:length(Stimulus(freqindex).RAttenIndices{attenindex})
				spikecol = Stimulus(freqindex).RAttenIndices{attenindex}(sweepindex);
				tmpcell{sweepindex} = Stimulus(freqindex).Spiketimes{unitindex, spikecol};
				%%% should probably subtract off Spikestart(sweepindex) to get
				%%% times relative to start of sweep %%%%%%%%%
				if ~isempty(tmpcell{sweepindex})
					tmpcell{sweepindex} = tmpcell{sweepindex} -  Stimulus(freqindex).Sweepstart(spikecol);
					tmpcell{sweepindex} = 0.001 * tmpcell{sweepindex};
				end
			end
			Spikes{unitindex, freqindex, attenindex} = tmpcell;
		end	% end of ATTENINDEX
	end	% end of FREQINDEX
end	% end of UNITINDEX

return

% plot a raster for each unit

unit = 1;
H = figure(1);

for freq = 1:Nstimuli
	
	for atten = 1:length(Stimulus(freqindex).RAttenVals)
		clear spiketimes_ms;
		for n = 1:length(Spikes{unit, freq, atten})
			spiketimes_ms{n} = 0.001 * Spikes{unit, freq, atten}(n);
		end
		rasterplot(spiketimes_ms, [0 1000], H);
		
		pause
	end
end

		

