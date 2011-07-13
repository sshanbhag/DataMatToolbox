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

% pre-allocate plot options
plotopts = struct( ...
	'time_limits',		[0 1000]		, ...
	'horizgap',			0.05			, ...
	'vertgap',			0.055			, ...
	'plotgap',			0.0125		, ...
	'filelabel',		D.Info.file	 ...
);

% create list of column labels (stimulus parameter - e.g., tone freq, wav name)
plotopts.columnlabels = cell(Nstimuli, 1);
for s = 1:Nstimuli
	if strcmp(StimList(s).Type, 'TONE')
		plotopts.columnlabels{s} = sprintf('Stim = %.2f', StimList(s).Var(1).values(1));
	elseif strcmp(StimList(col).Type, 'WAVFILE')
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

% for unit = 1:D.Info.Nunits
for unit = 1:1
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
	[Hplots, plotopts_out] = rasterpsthmatrix(tmpspikes, Natten, Nstimuli, plotopts);
end

%{
time_limits = [0 1000];

Natten = length(StimList(1).RAttenVals);

nrows = Natten;
ncols = Nstimuli;

horizgap = 0.05;
vertgap = 0.055;
plotgap = 0.0125;

plotwidth = (1 - ((ncols+1) * horizgap)) / ncols;
plotheight = (1 - ((nrows+2) * vertgap)) / (2*nrows);

pos1 = cell(nrows, ncols);
pos2 = cell(nrows, ncols);

for r = 1:nrows
	ypos1(r) = 1 - r*plotheight - (r-1)*plotheight - r*vertgap - (r-1)*plotgap;
	ypos2(r) = ypos1(r) - plotheight - plotgap;
	for c = 1:ncols
		xpos(c) = horizgap + ((c-1) * (plotwidth + horizgap));
		pos1{r, c} = [xpos(c) ypos1(r) plotwidth plotheight];
		pos2{r, c} = [xpos(c) ypos2(r) plotwidth plotheight];
	end
end

handles1 = cell(D.Info.Nunits, nrows, ncols);
handles2 = cell(D.Info.Nunits, nrows, ncols);

% for unit = 1:D.Info.Nunits
for unit = 1:1
	% for each unit, create a new page
	H = figure(unit);

	for atten = 1:Natten
		for var = 1:Nstimuli
			% select subplot location for rasters (pos1)
			subplot('Position', pos1{atten, var});
			% plot raster
			% store the axes handle returned by rasterplot in the handles2 cell array
			handles1{unit, atten, var} = rasterplot(Spikes{unit, var, atten}, time_limits);
			% turn off xtick labels, and turn off yticks
			set(gca, 'XTickLabel', []);
			set(gca, 'ytick', []);
			
			if (var == 1) & (atten == 1)
				unitstr = sprintf('Unit %d: ', unit);
			else
				unitstr = '';
			end

			if (atten == 1)
				if strcmp(StimList(var).Type, 'TONE')
					stimstr = sprintf('Stim = %.2f', StimList(var).Var(1).values(1));
				elseif strcmp(StimList(var).Type, 'WAVFILE')
					tmpstr = textscan(StimList(var).Var(2).values{1}, '%s', 'Delimiter', '\\');
					stimstr = ['Stim = ' tmpstr{1}{end}];
				else
					stimstr = ['Stim = ' num2str(StimList(var).Var(1).values(1))];
				end
			else
				stimstr = '';
			end
			
			if var == 1
				attenstr = sprintf('%d dB', StimList(var).RAttenVals(atten));
			else 
				attenstr = '';
			end
				
			titlestr = [unitstr ' ' stimstr ];
			title(titlestr, 'Interpreter', 'none')
			ylabel(attenstr, 'Interpreter', 'none')
	
			% select subplot location for psth (pos2)
			subplot('Position', pos2{atten, var});
			% build psth from spike data and plot using bar() function
			[histvals, bins] = psth(Spikes{unit, var, atten}, 5, 1000);
			bar(bins, histvals);
			% store the axes handle for the plot in the handles2 cell array
			handles2{unit, atten, var} = gca;
			% update time limits to match raster
			xlim(time_limits)
			% turn off x tick labels for all but the bottom row and
			% turn off y tick labels for all but the left column
			if atten ~= Natten
				set(gca, 'XTickLabel', []);
			end
			if var ~= 1
				set(gca, 'ytick', []);
			end
			
			% label the x-axis 'msec' on the lower left psth plot
			if (var == 1) && (atten == Natten)
				xlabel('msec')
			end
			% label the lower right plot axes with the input data file 
			if (var == Nstimuli) && (atten == Natten)
				xlabel(D.Info.file, 'Interpreter', 'none');
			end
				
		end
	end
end

%}


