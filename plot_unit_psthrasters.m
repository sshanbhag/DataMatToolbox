%-------------------------------------------------------------------
% plot_unit_psthrasters
%-------------------------------------------------------------------
% raster and psth data Plots

%% -----------------------------------------------------------------
%-------------------------------------------------------------------
% Some Constants
%-------------------------------------------------------------------
%-------------------------------------------------------------------

%-------------------------------------------------------------------
% limit to nonzero clusters?  1 = yes, 0 = no
%-------------------------------------------------------------------
NONZERO_LIMIT = 1;

%-------------------------------------------------------------------
% psth binsize
%-------------------------------------------------------------------
PSTHBINSIZE_MS = 10;
%-------------------------------------------------------------------
% mnemonics for conditions
%-------------------------------------------------------------------
PRE = 1;
MILD= 2;
CAT = 3;
%-------------------------------------------------------------------
% colors for different conditions
%-------------------------------------------------------------------
tcolors{PRE}	= 'k';
tcolors{MILD}	= 'b';
tcolors{CAT}	= 'r';
%-------------------------------------------------------------------
% background data parameters
%-------------------------------------------------------------------
BGTIME_MS = 5000;
BGSWEEPTIME_MS = 1000;
Nbgsweeps = round(BGTIME_MS / BGSWEEPTIME_MS);
bgwin = cell(Nbgsweeps, 1);
for b = 1:Nbgsweeps
	bgwin{b} = [b-1 b] * BGSWEEPTIME_MS;
end

%% -----------------------------------------------------------------
%-------------------------------------------------------------------
% Load data if necessary
%-------------------------------------------------------------------
%-------------------------------------------------------------------
inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
if NONZERO_LIMIT
	outpath = [inpath '/plots/nonzero_units'];
else
	outpath = [inpath '/plots/all_units'];
end
if ~exist(outpath, 'dir')
	mkdir(outpath);
end

if ~exist('bbnData', 'var')
	load(fullfile(inpath, 'BBNrate.mat'), 'validBBNList', 'bbnData', ...
						'spikeCountWindow', 'bg_spikeCountWindow')
end
if ~exist('lfhData', 'var')
	load(fullfile(inpath, 'LFHrate.mat'),  'validLFHList', ...
						'lfhData', 'spikeCountWindow', 'bg_spikeCountWindow')
end

Nwin = length(spikeCountWindow);
Nbgwin = length(bg_spikeCountWindow);
Nconditions = 3;
ntrials = zeros(Nconditions, 1);
tind = cell(Nconditions, 1);

%% -----------------------------------------------------------------
%-------------------------------------------------------------------
% bbnData plots
%-------------------------------------------------------------------
%-------------------------------------------------------------------
Ndata = length(validBBNList);

%-------------------------------------------------------------------
% plot background data
%-------------------------------------------------------------------
for d = 1:Ndata
	validDataIndices = validBBNList{d, 1};
	validUnitIndices = validBBNList{d, 5};
	[NvalidUnits, tmp] = size(validUnitIndices);
	
	if NONZERO_LIMIT
		% find indices to non-zero clusterd (i.e., sorted clusters only)
		clustn = zeros(NvalidUnits, Nconditions);
		for u = 1:NvalidUnits
			for c = 1:Nconditions
				clustn(u, c) = bbnData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
			end
		end	
		nonZeroUnits = find(clustn(:, 1));
		Nunits = length(nonZeroUnits);
	else
		nonZeroUnits = 1:NvalidUnits;
		Nunits = NvalidUnits;
	end
	
	% loop through units
	for u = 1:Nunits
		% get the index value for this non-zero unit
		uIndex = nonZeroUnits(u);
		
		tind{1} = 1:Nbgsweeps;
		if Nconditions > 1
			tc = length(tind{1});
			for c = 2:Nconditions
				tind{c} = tc + (1:Nbgsweeps);
				tc = tc + length(tind{c});
			end
		end
		
		S = cell(1, tc);
		
		for c = 1:Nconditions
			D = bbnData{validDataIndices(c)};
			% columns in validUnitIndices correspond to conditions/files
			U = D.UnitData(validUnitIndices(uIndex, c));
			tmpspikes =  U.BG_timestamps;
					
			tmpbg = cell(Nbgsweeps, 1);
			for trial = 1:Nbgsweeps
				spikewin = bgwin{trial};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < tmpspikes) & ...
											(tmpspikes < spikewin(2));
				% store spiketimes
				if ~isempty(validSpikes)
					tmpbg{trial} = tmpspikes(validSpikes) - spikewin(1);
				else
					tmpbg{trial} = [];
				end
			end
			S(tind{c}) = tmpbg;
		end

		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot raster
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hbgraster = figure(1);
		% select proper subplot for this unit
		subplot(1, Nunits, u);

		[h, hr] = rasterplot(S, [0 BGSWEEPTIME_MS], '.', 20);
		
		% set the colors for the data points
		for c = 1:Nconditions
			for t = tind{c}
				set(hr{t}, 'Color', tcolors{c});
			end
		end
		% labels/titles, etc.
		if u == 1
			tstr = {	validBBNList{d, 2}{1}, ...
						'BBN - Background', ...
					sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			xlabel('Time (ms)')
			set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')));
		else
			tstr = {	sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			set(gca, 'XTickLabel', []);
			set(gca, 'YTickLabel', []);
		end
		title(tstr, 'Interpreter', 'none');
		
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot PSTH
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hbgpsth = figure(2);
		for c = 1:Nconditions
			% select proper subplot for this unit
			subplot(Nconditions, Nunits, (c-1)*Nunits + u);
			[pcounts, pbins] = psth(S(tind{c}), PSTHBINSIZE_MS, BGSWEEPTIME_MS);
			bar(pbins, pcounts, 1, tcolors{c}, 'EdgeColor', tcolors{c});
			Hpsth = gcf;
			xlim([0 BGSWEEPTIME_MS]);
			if u == 1
				ylabel(sprintf('Condition %d', c), 'Color', tcolors{c});
				if c == 1
					tstr = {	validBBNList{d, 2}{1}, 'BBN - Background', ...
								sprintf('p%d c%d u%d',	U.UnitInfo.probe, ...
																U.UnitInfo.cluster, ...
																U.UnitInfo.unit ...
											) ...
								};
					title(tstr, 'Interpreter', 'none');
				elseif c == 3
					xlabel('Time (ms)');
				end
			elseif c == 1
				tstr = {	sprintf('p%d c%d u%d', ...
											U.UnitInfo.probe, ...
											U.UnitInfo.cluster, ...
											U.UnitInfo.unit	) ...
							};
				title(tstr, 'Interpreter', 'none');
				set(gca, 'XTickLabel', []);
			else
				set(gca, 'XTickLabel', []);
			end
		end

	end	% END uNITS LOOP
	drawnow
	
	% save plot to file
	outfile = [validBBNList{d, 2}{1} '_BBN_BGraster.fig'];
	fprintf('Writing plots...\n')
	fprintf('\traster plot to file: %s   ... ', outfile);
	saveas(Hbgraster, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');
	outfile = [validBBNList{d, 2}{1} '_BBN_BGpsth.fig'];
	fprintf('\tpsth plot to file: %s   ... ', outfile);
	saveas(Hbgpsth, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');	

end	% END dATA LOOP

%-------------------------------------------------------------------
% plot rasters for BBN sweeps
%-------------------------------------------------------------------
for d = 1:Ndata
	validDataIndices = validBBNList{d, 1};
	validUnitIndices = validBBNList{d, 5};
	validAttenIndices = validBBNList{d, 7};
	[NvalidUnits, tmp] = size(validUnitIndices);
	
	if NONZERO_LIMIT
		% find indices to non-zero clusterd (i.e., sorted clusters only)
		clustn = zeros(NvalidUnits, Nconditions);
		for u = 1:NvalidUnits
			for c = 1:Nconditions
				clustn(u, c) = bbnData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
			end
		end	
		nonZeroUnits = find(clustn(:, 1));
		Nunits = length(nonZeroUnits);
	else
		nonZeroUnits = 1:NvalidUnits;
		Nunits = NvalidUnits;
	end

	
	% loop through all units
	for u = 1:Nunits
		% get index to non-zero unit
		uIndex = nonZeroUnits(u);
		
		% get # of trials
		for c = 1:Nconditions
			attIndex = validAttenIndices(1, c);
			ntrials(c) = bbnData{validDataIndices(c)}.UnitData(validUnitIndices(uIndex, c)).ntrials(attIndex);
		end
		
		tind{1} = 1:ntrials(1);
		if Nconditions > 1
			tc = length(tind{1});
			for c = 2:Nconditions
				tind{c} = tc + (1:ntrials(c));
				tc = tc + length(tind{c});
			end
		end
		
		S = cell(1, tc);		
		for c = 1:Nconditions
			D = bbnData{validDataIndices(c)};
			attIndex = validAttenIndices(1, c);
			% columns in validUnitIndices correspond to conditions/files
			U = D.UnitData(validUnitIndices(uIndex, c));
			S(tind{c}) = U.spikes{attIndex};
		end
		
		sweepend =  fix(mean(1000*round(1e-6* diff(bbnData{validDataIndices(c)}.Stimulus.Sweepstart))));

		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot raster
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hraster = figure(3);
		% select proper subplot for this unit
		subplot(1, Nunits, u);
		
		[h, hr] = rasterplot(S, [0 sweepend], '.', 20);
		
		for c = 1:Nconditions
			for t = tind{c}
				set(hr{t}, 'Color', tcolors{c});
			end
		end
		
		if u == 1
			tstr = {	validBBNList{d, 2}{1}, ...
						'BBN', ...
					sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			xlabel('Time (ms)')
			set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')));
		else
			tstr = {	sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			set(gca, 'XTickLabel', []);
			set(gca, 'YTickLabel', []);
		end
		title(tstr, 'Interpreter', 'none')

		
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot PSTH
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hpsth = figure(4);
		for c = 1:Nconditions
			% select proper subplot for this unit
			subplot(Nconditions, Nunits, (c-1)*Nunits + u);
			[pcounts, pbins] = psth(S(tind{c}), PSTHBINSIZE_MS, BGSWEEPTIME_MS);
			bar(pbins, pcounts, 1, tcolors{c}, 'EdgeColor', tcolors{c});
			Hpsth = gcf;
			xlim([0 BGSWEEPTIME_MS]);
			if u == 1
				ylabel(sprintf('Condition %d', c), 'Color', tcolors{c});
				if c == 1
					tstr = {	validBBNList{d, 2}{1}, 'BBN', ...
								sprintf('p%d c%d u%d',	U.UnitInfo.probe, ...
																U.UnitInfo.cluster, ...
																U.UnitInfo.unit ...
											) ...
								};
					title(tstr, 'Interpreter', 'none');
				elseif c == 3
					xlabel('Time (ms)');
				end
			elseif c == 1
				tstr = {	sprintf('p%d c%d u%d', ...
											U.UnitInfo.probe, ...
											U.UnitInfo.cluster, ...
											U.UnitInfo.unit	) ...
							};
				title(tstr, 'Interpreter', 'none');
				set(gca, 'XTickLabel', []);
			else
				set(gca, 'XTickLabel', []);
			end
		end
		
		
	end	% END uNITS LOOP
	drawnow
	
	% save plot to file
	outfile = [validBBNList{d, 2}{1} '_BBN_raster.fig'];
	fprintf('Writing plots...\n')
	fprintf('\traster plot to file: %s   ... ', outfile);
	saveas(Hraster, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');
	outfile = [validBBNList{d, 2}{1} '_BBN_psth.fig'];
	fprintf('\tpsth plot to file: %s   ... ', outfile);
	saveas(Hpsth, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');	
	
end	% END dATA LOOP


%-------------------------------------------------------------------
%-------------------------------------------------------------------
%% lfhData
%-------------------------------------------------------------------
%-------------------------------------------------------------------
Ndata = length(validLFHList);

%-------------------------------------------------------------------
% plot background rasters
%-------------------------------------------------------------------
for d = 1:Ndata
	validDataIndices = validLFHList{d, 1};
	validUnitIndices = validLFHList{d, 5};
	[NvalidUnits, tmp] = size(validUnitIndices);
	
	if NONZERO_LIMIT
		% find indices to non-zero clusterd (i.e., sorted clusters only)
		clustn = zeros(NvalidUnits, Nconditions);
		for u = 1:NvalidUnits
			for c = 1:Nconditions
				clustn(u, c) = lfhData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
			end
		end	
		nonZeroUnits = find(clustn(:, 1));
		Nunits = length(nonZeroUnits);
	else
		nonZeroUnits = 1:NvalidUnits;
		Nunits = NvalidUnits;
	end
	
	% loop through units
	for u = 1:Nunits
		% get the index value for this non-zero unit
		uIndex = nonZeroUnits(u);
		
		tind{1} = 1:Nbgsweeps;
		if Nconditions > 1
			tc = length(tind{1});
			for c = 2:Nconditions
				tind{c} = tc + (1:Nbgsweeps);
				tc = tc + length(tind{c});
			end
		end
		
		S = cell(1, tc);
		
		for c = 1:Nconditions
			D = lfhData{validDataIndices(c)};
			% columns in validUnitIndices correspond to conditions/files
			U = D.UnitData(validUnitIndices(uIndex, c));
			tmpspikes =  U.BG_timestamps;
					
			tmpbg = cell(Nbgsweeps, 1);
			for trial = 1:Nbgsweeps
				spikewin = bgwin{trial};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < tmpspikes) & ...
											(tmpspikes < spikewin(2));
				% store spiketimes
				if ~isempty(validSpikes)
					tmpbg{trial} = tmpspikes(validSpikes) - spikewin(1);
				else
					tmpbg{trial} = [];
				end
			end
			S(tind{c}) = tmpbg;
		end

		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot raster
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hbgraster = figure(1);
		% select proper subplot for this unit
		subplot(1, Nunits, u);

		[h, hr] = rasterplot(S, [0 BGSWEEPTIME_MS], '.', 20);
		
		% set the colors for the data points
		for c = 1:Nconditions
			for t = tind{c}
				set(hr{t}, 'Color', tcolors{c});
			end
		end
		% labels/titles, etc.
		if u == 1
			tstr = {	validLFHList{d, 2}{1}, ...
						'LFH - Background', ...
					sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			xlabel('Time (ms)')
			set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')));
		else
			tstr = {	sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			set(gca, 'XTickLabel', []);
			set(gca, 'YTickLabel', []);
		end
		title(tstr, 'Interpreter', 'none');
		
		
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot PSTH
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hbgpsth = figure(2);
		for c = 1:Nconditions
			% select proper subplot for this unit
			subplot(Nconditions, Nunits, (c-1)*Nunits + u);
			[pcounts, pbins] = psth(S(tind{c}), PSTHBINSIZE_MS, BGSWEEPTIME_MS);
			bar(pbins, pcounts, 1, tcolors{c}, 'EdgeColor', tcolors{c});
			Hpsth = gcf;
			xlim([0 BGSWEEPTIME_MS]);
			if u == 1
				ylabel(sprintf('Condition %d', c), 'Color', tcolors{c});
				if c == 1
					tstr = {	validLFHList{d, 2}{1}, 'LFH - Background', ...
								sprintf('p%d c%d u%d',	U.UnitInfo.probe, ...
																U.UnitInfo.cluster, ...
																U.UnitInfo.unit ...
											) ...
								};
					title(tstr, 'Interpreter', 'none');
				elseif c == 3
					xlabel('Time (ms)');
				end
			elseif c == 1
				tstr = {	sprintf('p%d c%d u%d', ...
											U.UnitInfo.probe, ...
											U.UnitInfo.cluster, ...
											U.UnitInfo.unit	) ...
							};
				title(tstr, 'Interpreter', 'none');
				set(gca, 'XTickLabel', []);
			else
				set(gca, 'XTickLabel', []);
			end
		end
		

	end	% END uNITS LOOP
	drawnow
	
	% save plot to file
	outfile = [validLFHList{d, 2}{1} '_LFH_BGraster.fig'];
	fprintf('Writing plots...\n')
	fprintf('\traster plot to file: %s   ... ', outfile);
	saveas(Hbgraster, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');
	outfile = [validLFHList{d, 2}{1} '_LFH_BGpsth.fig'];
	fprintf('\tpsth plot to file: %s   ... ', outfile);
	saveas(Hbgpsth, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');	

	
end	% END dATA LOOP

%-------------------------------------------------------------------
% plot rasters for LFH sweeps
%-------------------------------------------------------------------
for d = 1:Ndata
	validDataIndices = validLFHList{d, 1};
	validUnitIndices = validLFHList{d, 5};
	validAttenIndices = validLFHList{d, 7};
	[NvalidUnits, tmp] = size(validUnitIndices);

	if NONZERO_LIMIT
		% find indices to non-zero clusterd (i.e., sorted clusters only)
		clustn = zeros(NvalidUnits, Nconditions);
		for u = 1:NvalidUnits
			for c = 1:Nconditions
				clustn(u, c) = lfhData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
			end
		end	
		nonZeroUnits = find(clustn(:, 1));
		Nunits = length(nonZeroUnits);
	else
		nonZeroUnits = 1:NvalidUnits;
		Nunits = NvalidUnits;
	end

	for u = 1:Nunits
		
		subplot(1, Nunits, u)

		uIndex = nonZeroUnits(u);
		
		% get # of trials
		for c = 1:Nconditions
			attIndex = validAttenIndices(1, c);
			ntrials(c) = lfhData{validDataIndices(c)}.UnitData(validUnitIndices(uIndex, c)).ntrials(attIndex);
		end
		
		tind{1} = 1:ntrials(1);
		if Nconditions > 1
			tc = length(tind{1});
			for c = 2:Nconditions
				tind{c} = tc + (1:ntrials(c));
				tc = tc + length(tind{c});
			end
		end
		
		S = cell(1, tc);
		
		for c = 1:Nconditions
			D = lfhData{validDataIndices(c)};
			attIndex = validAttenIndices(1, c);
			% columns in validUnitIndices correspond to conditions/files
			U = D.UnitData(validUnitIndices(uIndex, c));
			S(tind{c}) = U.spikes{attIndex};
		end
		
		sweepend =  fix(mean(1000*round(1e-6* diff(lfhData{validDataIndices(c)}.Stimulus.Sweepstart))));

		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot raster
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hraster = figure(3);
		% select proper subplot for this unit
		subplot(1, Nunits, u);
		% plot raster
		[h, hr] = rasterplot(S, [0 sweepend], '.', 20);
		
		for c = 1:Nconditions
			for t = tind{c}
				set(hr{t}, 'Color', tcolors{c});
			end
		end
		
		if u == 1
			tstr = {	validLFHList{d, 2}{1}, ...
						'LFH', ...
					sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			xlabel('Time (ms)');
			set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')));
		else
			tstr = {	sprintf('p%d c%d u%d', U.UnitInfo.probe, U.UnitInfo.cluster, U.UnitInfo.unit) };
			set(gca, 'XTickLabel', []);
			set(gca, 'YTickLabel', []);
		end
		title(tstr, 'Interpreter', 'none')
		
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% plot PSTH
		%-------------------------------------------------------------------
		%-------------------------------------------------------------------
		% select proper figure 
		Hpsth = figure(4);
		for c = 1:Nconditions
			% select proper subplot for this unit
			subplot(Nconditions, Nunits, (c-1)*Nunits + u);
			[pcounts, pbins] = psth(S(tind{c}), PSTHBINSIZE_MS, BGSWEEPTIME_MS);
			bar(pbins, pcounts, 1, tcolors{c}, 'EdgeColor', tcolors{c});
			Hpsth = gcf;
			xlim([0 BGSWEEPTIME_MS]);
			if u == 1
				ylabel(sprintf('Condition %d', c), 'Color', tcolors{c});
				if c == 1
					tstr = {	validLFHList{d, 2}{1}, 'LFH', ...
								sprintf('p%d c%d u%d',	U.UnitInfo.probe, ...
																U.UnitInfo.cluster, ...
																U.UnitInfo.unit ...
											) ...
								};
					title(tstr, 'Interpreter', 'none');
				elseif c == 3
					xlabel('Time (ms)');
				end
			elseif c == 1
				tstr = {	sprintf('p%d c%d u%d', ...
											U.UnitInfo.probe, ...
											U.UnitInfo.cluster, ...
											U.UnitInfo.unit	) ...
							};
				title(tstr, 'Interpreter', 'none');
				set(gca, 'XTickLabel', []);
			else
				set(gca, 'XTickLabel', []);
			end
		end

	end	% END uNITS LOOP
	drawnow
	

	% save plot to file
	outfile = [validLFHList{d, 2}{1} '_LFH_raster.fig'];
	fprintf('Writing plots...\n')
	fprintf('\traster plot to file: %s   ... ', outfile);
	saveas(Hraster, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');
	outfile = [validLFHList{d, 2}{1} '_LFH_psth.fig'];
	fprintf('\tpsth plot to file: %s   ... ', outfile);
	saveas(Hpsth, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');	

end	% END dATA LOOP









