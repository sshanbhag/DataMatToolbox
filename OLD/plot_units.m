%% data Plots

Nwin = length(spikeCountWindow);
Nbgwin = length(bg_spikeCountWindow);
Nconditions = 3;
ntrials = zeros(Nconditions, 1);
tind = cell(Nconditions, 1);

PRE = 1;
MILD= 2;
CAT = 3;

tcolors = {'b', 'g', 'r'};

inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
outpath = [inpath '/rasters/nonzero_units'];
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

Ndata = length(validBBNList);

for d = 1:Ndata
	validDataIndices = validBBNList{d, 1};
	validUnitIndices = validBBNList{d, 5};
	validAttenIndices = validBBNList{d, 7};
	[NvalidUnits, tmp] = size(validUnitIndices);
	
	
	clustn = zeros(NvalidUnits, Nconditions);
	for u = 1:NvalidUnits
		for c = 1:Nconditions
			clustn(u, c) = bbnData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
		end
	end
	
	nonZeroUnits = find(clustn(:, 1));
	Nunits = length(nonZeroUnits);
	

	for u = 1:Nunits
		
		subplot(1, Nunits, u)
		
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

	end	% END uNITS LOOP
	drawnow
	
	% save plot to file
	outfile = [validBBNList{d, 2}{1} '_BBN_raster.fig'];
	fprintf('Writing plot to file: %s   ... ', outfile);
	saveas(gcf, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');
	
end	% END dATA LOOP


Ndata = length(validLFHList);
for d = 1:Ndata
	validDataIndices = validLFHList{d, 1};
	validUnitIndices = validLFHList{d, 5};
	validAttenIndices = validLFHList{d, 7};
	[NvalidUnits, tmp] = size(validUnitIndices);

	clustn = zeros(NvalidUnits, Nconditions);
	for u = 1:NvalidUnits
		for c = 1:Nconditions
			clustn(u, c) = lfhData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
		end
	end
	
	nonZeroUnits = find(clustn(:, 1));
	Nunits = length(nonZeroUnits);
	

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

	end	% END uNITS LOOP
	drawnow
	
	outfile = [validLFHList{d, 2}{1} '_LFH_raster.fig'];
	fprintf('Writing plot to file: %s   ... ', outfile);
	saveas(gcf, fullfile(outpath, outfile), 'fig');
	fprintf('... done.\n');
	
end	% END dATA LOOP









