function varargout = plotData(D, Stimulus)
%------------------------------------------------------------------------
% plotData
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% Input Arguments:
%	D		Datawave struct
% 	Stimulus	Stimulus struct
% 
% Output Arguments:
%		Hplots	figure handles
%
%------------------------------------------------------------------------
% See:  
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 	4 July 2011 (SJS):
%
% Revisions:
%	-	28 July, 2011 (SJS): functionalized
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;


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

if toneIndex
	StimList = ToneStim;
	Nstimuli = toneIndex;
elseif wavIndex
	StimList = WavStim;
	Nstimuli = wavIndex;
elseif noiseIndex
	StimList = NoiseStim;
	Nstimuli = noiseIndex;
elseif unknownIndex
	StimList = UnknownStim;
	Nstimuli = unknownIndex;
else
	error('%s: no known stimulus types found', mfilename);
end

% find ranges of attenuation
for s = 1:Nstimuli
	% R channel attenuation
	AttenR{s} = StimList(s).RAttenVals;
	AttenR_indices{s} = StimList(s).RAttenIndices;
end


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
	elseif strcmp(StimList(s).Type, 'NOISE')
		plotopts.columnlabels{s} = sprintf('Stim = BBN');
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
	Hfig(unit) = figure;
	
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


if nargout
	varargout{1} = Hplots;
end