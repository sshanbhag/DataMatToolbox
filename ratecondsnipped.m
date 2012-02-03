



















































%{
	%--------------------------------------------------------------------------
	% construct data arrays for analysis
	%--------------------------------------------------------------------------
	% use buildSpikes function to build the spikes arrays
	%--------------------------------------------------------------------------
	BBN.Spikes = buildSpikes(BBN.Stimulus);
	LFH.Spikes = buildSpikes(LFH.Stimulus);

	%--------------------------------------------------------------------------
	% determine number of trials (use second atten value to eliminate extra
	% trial at end)
	%--------------------------------------------------------------------------
	BBN.ntrials = length(BBN.Spikes{1, 1, minCommonAttenIndices(1)});
	LFH.ntrials = length(LFH.Spikes{1, 1, minCommonAttenIndices(2)});

	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% loop through units
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% clear UnitData variable
	clear UnitData;

	% initialize the unitIndex counter - this will be used to index the analyzed
	% data output in MUdata
	dIndex = 0;

	% loop through the units to compare for this pair of recordings (files)
	for u = 1:Nunits
		% increment counter
		dIndex = dIndex + 1;

		% get the 1 X 2 unitNumber information [BBNunit#, LFHunit#]
		unitNumber = UnitIndex(u, :);

		% store current unit information in UnitData structure
		UnitData(u).UnitNumber = unitNumber;
		UnitData(u).UnitInfo(1) = BBN.D.Info.UnitInfo(unitNumber(1));
		UnitData(u).UnitInfo(2) = LFH.D.Info.UnitInfo(unitNumber(2));
		UnitData(u).UnitList = UnitList(u, :);

		% get all spikes for current unit and attenuation, condition A
		UnitData(u).spikes{1} = BBN.Spikes{unitNumber(1), 1, minCommonAttenIndices(1)};
		UnitData(u).spikes{2} = LFH.Spikes{unitNumber(2), 1, minCommonAttenIndices(2)};

%}


fp = fopen(fullfile(outpath, outcsvfile), 'w');

fprintf(fp, 'SpikeCountWindowLength,%f,msec,\n', RESPONSEWINDOW_MS);

fprintf(fp, 'BBN_File,');
fprintf(fp, 'LFH_File,');
fprintf(fp, 'BBN_UnitNum,');
fprintf(fp, 'LFH_UnitNum,');
fprintf(fp, 'BBN_Probe,');
fprintf(fp, 'LFH_Probe,');
fprintf(fp, 'BBN_Cluster,');
fprintf(fp, 'LFH_Cluster,');
fprintf(fp, 'BBN_BGcount,');
fprintf(fp, 'LFH_BGcount,');
fprintf(fp, 'BBN_Atten_dB,');
fprintf(fp, 'LFH_Atten_dB,');
fprintf(fp, 'BBN_Net_mean,');
fprintf(fp, 'LFH_Net_mean,');
fprintf(fp, 'BBN_Net_sd,');
fprintf(fp, 'LFH_Net_sd,');
fprintf(fp, 'BBN_Count_mean,');
fprintf(fp, 'LFH_Count_mean,');
fprintf(fp, 'BBN_Count_sd,');
fprintf(fp, 'LFH_Count_sd,');
fprintf(fp, '\n');

rindx = 0;
for f = 1:length(validBBNList)
	Nunits = length(Data{f}.UnitData);
	tmpD = Data{f};
	for u = 1:Nunits
		tmpU = Data{f}.UnitData(u);
		fprintf(fp, '%s,', tmpD.BBNinfo.file);
		fprintf(fp, '%s,', tmpD.LFHinfo.file);
		fprintf(fp, '%d,', tmpU.UnitInfo(1).unit);
		fprintf(fp, '%d,', tmpU.UnitInfo(2).unit);
		fprintf(fp, '%d,', tmpU.UnitInfo(1).probe);
		fprintf(fp, '%d,', tmpU.UnitInfo(2).probe);
		fprintf(fp, '%d,', tmpU.UnitInfo(1).cluster);
		fprintf(fp, '%d,', tmpU.UnitInfo(2).cluster);
		fprintf(fp, '%.4f,', tmpU.BG_count(1));
		fprintf(fp, '%.4f,', tmpU.BG_count(2));
		fprintf(fp, '%.4f,', tmpD.minAtten(1));
		fprintf(fp, '%.4f,', tmpD.minAtten(2));
		fprintf(fp, '%.4f,', tmpU.Net_mean(1));
		fprintf(fp, '%.4f,', tmpU.Net_mean(2));
		fprintf(fp, '%.4f,', tmpU.Net_std(1));
		fprintf(fp, '%.4f,', tmpU.Net_std(2));
		fprintf(fp, '%.4f,', tmpU.Count_mean{1});
		fprintf(fp, '%.4f,', tmpU.Count_mean{2});
		fprintf(fp, '%.4f,', tmpU.Count_std{1});
		fprintf(fp, '%.4f,', tmpU.Count_std{2});
		fprintf(fp, '\n');
		
		rindx = rindx + 1;
		c = 1;
		Rdata{rindx, c} = tmpD.BBNinfo.file;	c = c + 1;
		Rdata{rindx, c} = tmpD.LFHinfo.file;	c = c + 1;
		Rdata{rindx, c} = tmpU.UnitInfo(1).unit;	c = c + 1;
		Rdata{rindx, c} = tmpU.UnitInfo(2).unit;	c = c + 1;
		Rdata{rindx, c} = tmpU.UnitInfo(1).probe;	c = c + 1;
		Rdata{rindx, c} = tmpU.UnitInfo(2).probe;	c = c + 1;
		Rdata{rindx, c} = tmpU.UnitInfo(1).cluster;	c = c + 1;
		Rdata{rindx, c} = tmpU.UnitInfo(2).cluster;	c = c + 1;
		Rdata{rindx, c} = tmpU.BG_count(1);	c = c + 1;
		Rdata{rindx, c} = tmpU.BG_count(2);	c = c + 1;
		Rdata{rindx, c} = tmpD.minAtten(1);	c = c + 1;
		Rdata{rindx, c} = tmpD.minAtten(2);	c = c + 1;
		Rdata{rindx, c} = tmpU.Net_mean(1);	c = c + 1;
		Rdata{rindx, c} = tmpU.Net_mean(2);	c = c + 1;
		Rdata{rindx, c} = tmpU.Net_std(1);	c = c + 1;
		Rdata{rindx, c} = tmpU.Net_std(2);	c = c + 1;
		Rdata{rindx, c} = tmpU.Count_mean{1};	c = c + 1;
		Rdata{rindx, c} = tmpU.Count_mean{2};	c = c + 1;
		Rdata{rindx, c} = tmpU.Count_std{1};	c = c + 1;
		Rdata{rindx, c} = tmpU.Count_std{2};	c = c + 1;
	end
end

if fp ~= 1
	fclose(fp);
end


%--------------------------------------------------------------------------
% write to mat file
%--------------------------------------------------------------------------
clear BBN
clear LFH
BBN.file =			Rdata(:, 1);
LFH.file =			Rdata(:, 2);
BBN.unit =			cell2mat(Rdata(:, 3));
LFH.unit =			cell2mat(Rdata(:, 4));
BBN.probe =			cell2mat(Rdata(:, 5));
LFH.probe =			cell2mat(Rdata(:, 6));
BBN.cluster =		cell2mat(Rdata(:, 7));
LFH.cluster =		cell2mat(Rdata(:, 8));
BBN.BG_count =		cell2mat(Rdata(:, 9));
LFH.BG_count =		cell2mat(Rdata(:, 10));
BBN.minAtten =		cell2mat(Rdata(:, 11));
LFH.minAtten =		cell2mat(Rdata(:, 12));
BBN.Net_mean =		cell2mat(Rdata(:, 13));
LFH.Net_mean =		cell2mat(Rdata(:, 14));
BBN.Net_sd =		cell2mat(Rdata(:, 15));
LFH.Net_sd =		cell2mat(Rdata(:, 16));
BBN.Count_mean =	cell2mat(Rdata(:, 17));
LFH.Count_mean =	cell2mat(Rdata(:, 18));
BBN.Count_sd =		cell2mat(Rdata(:, 19));
LFH.Count_sd =		cell2mat(Rdata(:, 20));

BBN.BGrate = BBN.BG_count ./ 5;
LFH.BGrate = LFH.BG_count ./ 5;
BBN.Netrate = BBN.Net_mean ./ 1;
LFH.Netrate = LFH.Net_mean ./ 1;
BBN.Rate_mean = BBN.Count_mean ./ (0.001 * RESPONSEWINDOW_MS);
LFH.Rate_mean = LFH.Count_mean ./ (0.001 * RESPONSEWINDOW_MS);

save(fullfile(outpath, outcountmatfile), 'Rdata', 'BBN', 'LFH', 'Data', '-MAT')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plot
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% find non-zero clusters
nonz = find(BBN.cluster > 0);

% plot counts information
figure
subplot(131)
plot(BBN.Net_mean(nonz), LFH.Net_mean(nonz), '.')
title('Net # spikes')
xlabel('BBN');
ylabel('LFH');
xm = xlim;
ym = ylim;
maxval = max([xm(2) ym(2)]);
xlim([0 maxval])
ylim([0 maxval])
axis square
line(xlim, ylim)

subplot(132)
plot(BBN.Count_mean(nonz), LFH.Count_mean(nonz), '.')
title(sprintf('Window Mean (%d ms) # spikes', RESPONSEWINDOW_MS))
xlabel('BBN');
ylabel('LFH');
xm = xlim;
ym = ylim;
maxval = max([xm(2) ym(2)]);
xlim([0 maxval])
ylim([0 maxval])
axis square
line(xlim, ylim)

subplot(133)
plot(BBN.BG_count(nonz), LFH.BG_count(nonz), '.')
title('Background # spikes')
xlabel('BBN');
ylabel('LFH');
xm = xlim;
ym = ylim;
maxval = max([xm(2) ym(2)]);
xlim([0 maxval])
ylim([0 maxval])
axis square
line(xlim, ylim)


figure
subplot(131)
plot(BBN.Netrate(nonz), LFH.Netrate(nonz), '.')
title('Net rate')
xlabel('BBN');
ylabel('LFH');
xm = xlim;
ym = ylim;
maxval = max([xm(2) ym(2)]);
xlim([0 maxval])
ylim([0 maxval])
axis square
line(xlim, ylim)

subplot(132)
plot(BBN.Rate_mean(nonz), LFH.Rate_mean(nonz), '.')
title(sprintf('Window Mean (%d ms) rate', RESPONSEWINDOW_MS))
xlabel('BBN');
ylabel('LFH');
xm = xlim;
ym = ylim;
maxval = max([xm(2) ym(2)]);
xlim([0 maxval])
ylim([0 maxval])
axis square
line(xlim, ylim)

subplot(133)
plot(BBN.BGrate(nonz), LFH.BGrate(nonz), '.')
title('Background rate')
xlabel('BBN');
ylabel('LFH');
xm = xlim;
ym = ylim;
maxval = max([xm(2) ym(2)]);
xlim([0 maxval])
ylim([0 maxval])
axis square
line(xlim, ylim)


% plot histograms

figure
subplot(411)
hist(BBN.Net_mean(nonz), 30);
title('BBN Net Mean # spikes')
subplot(412)
hist(LFH.Net_mean(nonz), 30)
title('LFH Net Mean # spikes')
subplot(413)
hist(BBN.Count_mean(nonz), 30); 
title(sprintf('BBN Window Mean (%d ms) # spikes', RESPONSEWINDOW_MS))
subplot(414)
hist(LFH.Count_mean(nonz), 30)
title(sprintf('LFH Window Mean (%d ms) # spikes', RESPONSEWINDOW_MS))


