%------------------------------------------------------------------------
% plotFRAData
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
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 	19 July 2011 (SJS):
%	-	split off from plotData to handle the FRA (frequency response area)
%		data that are rather different from the BBN (broad band (sometimes)
% 		noise) and WAV (.wav file) stimuli.
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
			ToneFreq{toneIndex} = Stimulus(s).Var(1).values;
		case 'WAVFILE'
			disp(['Stimulus #' num2str(s) ' is WAVFILE, skipping...']);
		case 'NOISE'
			disp(['Stimulus #' num2str(s) ' is NOISE, skipping...']);
		otherwise
			disp(['Stimulus #' num2str(s) ' is UNKNOWN, skipping...']);
	end
end

if toneIndex
	StimList = ToneStim;
	Nstimuli = toneIndex;
else
	error('%s: no TONE stimulus types found in data!', mfilename);
end

% find frequencies and sort stimuli
clear tmp;
for s = 1:Nstimuli
	for n = 1:length(StimList(s).Var)
 		if strncmp(StimList(s).Var(n).name, 'ToneFreqR', length('ToneFreqR'))
 			tmp(s) = StimList(s).Var(n).values(1);
		end
	end
end

% count number of different tone frequencies
Nfreqs = length(tmp);
% sort freqs from low to high, keeping indices
[Freqs, findex] = sort(tmp);
% use indices to sort stimuli in StimList from low to high frequency
StimList = StimList(findex);

% find ranges of attenuation
for s = 1:Nstimuli
	% R channel attenuation
	AttenR{s} = StimList(s).RAttenVals;
	AttenR_indices{s} = StimList(s).RAttenIndices;
end

% count number of different attenuation settings (assume same for all
% stimuli)
Natten = length(StimList(1).RAttenVals);
% assume that all stimuli are presented at same intensity levels -
% should probably do a check on this
AttenVals = StimList(1).RAttenVals;


%----------------------------------------------------------------------
%----------------------------------------------------------------------
% build Spikes cell array to hold spike time information for 
%----------------------------------------------------------------------
% Spikes has three dimensions: 
%		(1) unit index
%		(2) stimulus index (varying parameter, e.g., wav file, tone freq)
%		(3) attenuation value
%----------------------------------------------------------------------
%----------------------------------------------------------------------
Spikes = buildSpikes(StimList);

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% plot!
%----------------------------------------------------------------------
%----------------------------------------------------------------------
FreqAttenSpikeCount = cell(D.Info.Nunits, 1);

for unit = 1:D.Info.Nunits
	SpikeCount = zeros(Natten, Nfreqs);

	% loop through frequencies
	for freq = 1:Nfreqs	
		% loop through attenuation from low to high
		fatten = 0;
		for atten = 1:Natten
			fatten = fatten + 1;
			tmp = Spikes{unit, freq, atten};

			% add up number of spikes for this freq and atten combination
			for rep = 1:length(tmp)
				if ~isempty(tmp{rep})
					SpikeCount(fatten, freq) = SpikeCount(atten, freq) + length(tmp{rep});
				end
			end
			
		end
	end

	figure
	subplot(211)
	pcolor(log10(Freqs), AttenVals, SpikeCount)
	colorbar
	xlabel('Log Frequency (kHz)')
	ylabel('Attenuation (dB)')
	title(sprintf('Unit %d: Spike Count', unit));
	xt = get(gca, 'XTick')
	xtl = cell(length(xt), 1);
	get(gca, 'XTickLabel')
	for n = 1:length(xt)
		xtl{n} = sprintf('%.0f', 0.001 * 10^xt(n));
	end
	set(gca, 'XTickLabel', xtl);
	
	
	subplot(212)
	waterfall(log10(Freqs), AttenVals, SpikeCount);
	xlabel('Log Frequency (kHz)')
	ylabel('Attenuation (dB)')
	zlabel('Spike Count');
	xt = get(gca, 'XTick')
	xtl = cell(length(xt), 1);
	get(gca, 'XTickLabel')
	for n = 1:length(xt)
		xtl{n} = sprintf('%.0f', 0.001 * 10^xt(n));
	end
	set(gca, 'XTickLabel', xtl);

	FreqAttenSpikeCount{unit} = SpikeCount;
	
end






%{

%----------------------------------------------------------------------
%----------------------------------------------------------------------
% need to get max inter stimulus interval for determining plot limits

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
	[Hplots(unit), Plotopts(unit)] = rasterpsthmatrix(tmpspikes, plotopts);
end


%}