%------------------------------------------------------------------------
% spikerater.m
%------------------------------------------------------------------------
% script for spike rate computation and analysis from DataMat toolbox
%------------------------------------------------------------------------
% Some Notes:
%------------------------------------------------------------------------

clear all

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some constants/settings
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% size of time window after stimulus onset in which to count spikes
RESPONSEWINDOW_MS = 250;

% default sound onset time time
soundOnsetTime = 0;
spikeCountWindow{1} = [0 RESPONSEWINDOW_MS];
Nwin = length(spikeCountWindow);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some options
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles';

types_to_process = {'BBN', 'LFH'};

outpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles';
matrootname = 'Data_BBN_LFH'
outmatfile = [matrootname '.mat'];
outcountmatfile = sprintf('%s_%dms.mat', matrootname, RESPONSEWINDOW_MS);
outcsvfile = sprintf('%s_%dms.csv', matrootname, RESPONSEWINDOW_MS);



%--------------------------------------------------------------------------
% define some plot constants/settings
%--------------------------------------------------------------------------


if exist(fullfile(outpath, outmatfile), 'file')
	load(fullfile(outpath, outmatfile));
else

	%------------------------------------------------------------------------
	% List of files and units to analyze
	%------------------------------------------------------------------------

	% 1) need to find matching pairs of files

	% get all BBN files first, then search for matching LFH files?

	% get list of files, store info in fstruct struct array
	fstruct = dir(inpath);

	% make sure there are files
	if isempty(fstruct)
		error('%s: no files found in %s', mfilename, inpath);
	end

	% find files that aren't empty
	tmp = struct2cell(fstruct);
	bytes = cell2mat(tmp(3, :));
	clear tmp

	% store valid files
	[tmp, valid_findex] =  find(bytes > 0);
	nfiles = length(valid_findex);

	bbnCount = 0;
	lfhCount = 0;
	% loop through valid files
	for findex = 1:nfiles
		% build full filename (with path)
		fname = fstruct(valid_findex(findex)).name;

		% check if this is one of the types to process
		if strfind(fname, 'BBN')
			bbnCount = bbnCount + 1;
			bbnFiles{bbnCount} = fname;
		elseif strfind(fname, 'LFH')
			lfhCount = lfhCount + 1;
			lfhFiles{lfhCount} = fname;
		end
	end

	for n = 1:bbnCount
		bbn_idstr{n} = regexprep(bbnFiles{n}, ('_BBN'), '');
	end

	for n = 1:lfhCount
		lfh_idstr{n} = regexprep(lfhFiles{n}, ('_LFH'), '');
	end

	match_index = 0;
	for n = 1:bbnCount
		lfhsrch = strcmp(bbn_idstr{n}, lfh_idstr);

		if sum(lfhsrch)
			match_index = match_index + 1;
			bbnlfhMatch{match_index, 1} = n;		
			bbnlfhMatch{match_index, 2} = find(lfhsrch);
		end

	end

	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% Loop through file pairs
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% Two issues needed to be dealt with:
	% (1) # of units in each file (from each "test" with a specific stimulus, 
	% 	e.g., BBN and LFH) is not guaranteed to be the same.  Thus, the 
	% 	unit index from loadDWfile output will not necessarily match.
	% 	
	% 	However, each unit's probe # and cluster # will be consistent between
	% 	files.  This unit information is stored in the UnitList structure and can
	% 	be compared between the two files and the common units identified for 
	% 	further analysis.
	% 	
	% (2) The attenuation values used for the stimuli in each file may, or may not
	% 	 match.
	% 	 
	% 	As a workaround, the lowest common attenuation value is found and analysis 
	% 	will be done on those data
	% 
	%--------------------------------------------------------------------------

	% pre-allocate Data cell array to hold data for each pair of BBN and LFH files
	Data = cell(1, match_index);

	% loop through file list
	for fIndx = 1:match_index
		%------------------------------------------------------------------------
		% Load Data files
		%------------------------------------------------------------------------
		bbnfile = bbnFiles{bbnlfhMatch{fIndx, 1}};
		lfhfile = lfhFiles{bbnlfhMatch{fIndx, 2}};
		BBN = load(fullfile(inpath, bbnfile), 'D', 'Stimulus', 'Background');
		LFH = load(fullfile(inpath, lfhfile), 'D', 'Stimulus', 'Background');

		fprintf('%d\n', fIndx);
		fprintf('\t%s\n',	BBN.D.Info.file);
		fprintf('\t%s\n\n',	LFH.D.Info.file);

		%------------------------------------------------------------------------
		% make sure # of units match
		%------------------------------------------------------------------------
		%	store matching units information in UnitList
		%------------------------------------------------------------------------
		if BBN.D.Info.Nunits ~= LFH.D.Info.Nunits
			warning('%s: # Units mismatch!!!!!!!!!!!\n', mfilename)
			fprintf('\tBBN file has %d\tunits\n', BBN.D.Info.Nunits);
			fprintf('\tLFH file has %d\tunits\n\n', LFH.D.Info.Nunits);

			% make some local copies of the UnitList struct from each file
			blist = BBN.D.Info.UnitList;
			fprintf('\tBBN Unit Info:\n')
			fprintf('\t\tUnit\t\tProbe\t\tCluster\n')
			for p = 1:BBN.D.Info.Nunits
				fprintf('\t\t%d\t\t\t%d\t\t\t%d\n', blist(p, 3), blist(p, 1), blist(p, 2) );
			end

			llist = LFH.D.Info.UnitList;
			fprintf('\tLFH Unit Info:\n')
			fprintf('\t\tUnit\t\tProbe\t\tCluster\n')
			for p = 1:LFH.D.Info.Nunits
				fprintf('\t\t%d\t\t\t%d\t\t\t%d\n', llist(p, 3), llist(p, 1), llist(p, 2) );
			end

			% if mismatch # of units, find units that are common to both
			% i.e., same probe and cluster #
			% store this list in UnitList which is a n X 3 matrix with
			%  [ Probe#		Cluster#	] as columns
			% and store unit numbers in 
			% UnitIndex[BBN.D.Info.UnitList(:, 3) LFH.D.Info.UnitList(:, 3)]
			if BBN.D.Info.Nunits > LFH.D.Info.Nunits
				% total # of matching units MUST be the lower of the 
				% # of units in the two files... in this case, LFH
				UnitList = zeros(LFH.D.Info.Nunits, 2);
				UnitIndex = UnitList;
				u = 0;
				for p = 1:BBN.D.Info.Nunits
					lFlag = 0;
					for l = 1:LFH.D.Info.Nunits
						% check if probe # (col 1) and cluster # (col 2) match
						if all(blist(p, 1:2) == llist(l, 1:2))
							lFlag = l;
						end
					end
					% if match found, store the unit information
					if lFlag
						u = u + 1;
						UnitList(u, :) = blist(p, 1:2);
						UnitIndex(u, :) = [p lFlag];
					end
				end
			else
				% fewer units in BBN file, so that is the # of matching units
				UnitList = zeros(BBN.D.Info.Nunits, 2);
				UnitIndex = UnitList;
				u = 0;
				for p = 1:LFH.D.Info.Nunits
					bFlag = 0;
					for l = 1:BBN.D.Info.Nunits
						if all(llist(p, 1:2) == blist(l, 1:2))
							bFlag = l;
						end
					end

					if bFlag
						u = u + 1;
						UnitList(u, :) = llist(p, 1:2);
						UnitIndex(u, :) = [bFlag p];
					end
				end
			end

			Nunits = length(UnitList(:, 1));
		else
			% units match so set UnitList to one of the lists from the 
			% D structures (for no particular reason, BBN is used)
			Nunits = BBN.D.Info.Nunits;
			UnitList = BBN.D.Info.UnitList(:, 1:2);
			UnitIndex = [BBN.D.Info.UnitList(:, 3) LFH.D.Info.UnitList(:, 3)];
		end

		%--------------------------------------------------------------------------
		% determine number of Attenuation levels by looking at the Stimulus.Var
		% information in each condition's data structure
		%--------------------------------------------------------------------------
		% check number of attenuation values

		% # of attenuation values in BBN file
		nVarA = length(BBN.Stimulus.Var);
		attenindex = 0;
		for n = 1:nVarA
			% find R channel attenuation
			if strcmp(BBN.Stimulus.Var(n).name, 'AttenuationR')
				attenindex = n;
			end
		end
		% if found, store the values
		if attenindex
			attenValsA = BBN.Stimulus.Var(attenindex).values;
		else
			% if not, error
			error('%s: attenuation not found in Var list for file %s', mfilename, bbnfile);
		end
		NattenA = length(attenValsA);

		% repeat process for LFH file
		nVarB = length(LFH.Stimulus.Var);
		attenindex = 0;
		for n = 1:nVarB
			if strcmp(LFH.Stimulus.Var(n).name, 'AttenuationR')
				attenindex = n;
			end
		end
		if attenindex
			attenValsB = LFH.Stimulus.Var(attenindex).values;
		else
			error('%s: attenuation not found in Var list for file %s', mfilename, lfhfile);
		end
		NattenB = length(attenValsB);

		% find minimum common attenuation value
		ccount = 0;
		commonAtten = [];
		% loop through attenuation values
		for n = 1:NattenA
			% find matches
			cind = find(attenValsA(n) == attenValsB);
			if ~isempty(cind)
				% store matching values
				ccount = ccount + 1;
				commonAtten(ccount) = attenValsA(n);
			end
		end
		% if found, find the minimum common value
		if ccount
			minCommonAtten = min(commonAtten);
			minCommonAttenIndices = ...
				[find(attenValsA == minCommonAtten) find(attenValsB == minCommonAtten)];
		else
			% if none found, use the lowest one from each file
			warning('%s: no common attenuation values in files %s and %s', mfilename, bbnfile, lfhfile);
			fprintf('\tUsing lowest value from each file\n')
			minCommonAtten = [min(attenValsA) min(attenValsB)];
			minCommonAttenIndices = ...
				[find(attenValsA == minCommonAtten(1)) find(attenValsB == minCommonAtten(2))];
			fprintf('\tUsing lowest value from each file\n')
			fprintf('\t\tattenValA = %d\n', minCommonAtten(1));
			fprintf('\t\tattenValB = %d\n', minCommonAtten(2));
		end

		%------------------------------------------------------------------------
		% List of units to analyze
		%------------------------------------------------------------------------
		fprintf('\nCommon units found: %d\n', Nunits);
		fprintf('\t\tProbe\t\tCluster\n')
		for p = 1:Nunits
			fprintf('\t\t\t%d\t\t\t%d\n', UnitList(p, 1), UnitList(p, 2));
		end
		fprintf('Using Attenuation %d dB\n', minCommonAtten);

		%--------------------------------------------------------------------------
		% now, need to parse filenames in order to extract condition
		% 		...BBN_..._1.mat  ---> pre odor (clean bedding)
		% 		...BBN_..._2.mat  ---> mild aversive odor (clean cotton)
		% 		...BBN_..._3.mat  ---> strong aversive odof (cat odor)
		%--------------------------------------------------------------------------

		for f = 1:length(Data)
			

		end



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

			% Assign spike counts to the (surprise) SpikeCounts array
			% assuming that both conditions have same number of trials... 
			% pre allocate SpikeCounts matrix
			NetSpikeCounts{1} = zeros(1, BBN.ntrials);
			NetSpikeCounts{2} = zeros(1, LFH.ntrials);
			for trial = 1:BBN.ntrials
				NetSpikeCounts{1}(trial) = length(UnitData(u).spikes{1}{trial});
			end
			for trial = 1:LFH.ntrials
				NetSpikeCounts{2}(trial) = length(UnitData(u).spikes{2}{trial});
			end

			% store in overall data structure
			UnitData(u).NetCounts = NetSpikeCounts;

			%------------------------------------------------------------------------
			% get background spikes
			%------------------------------------------------------------------------
			UnitData(u).BG_timestamps{1} = BBN.Background.Spiketimes{unitNumber(1)};
			UnitData(u).BG_count(1) = length(BBN.Background.Spiketimes{unitNumber(1)});
			UnitData(u).BG_timestamps{2} = LFH.Background.Spiketimes{unitNumber(2)};
			UnitData(u).BG_count(2) = length(LFH.Background.Spiketimes{unitNumber(2)});

		end	% END OF u LOOP

		%--------------------------------------------------------------------------
		%--------------------------------------------------------------------------
		% process data
		%--------------------------------------------------------------------------
		%--------------------------------------------------------------------------
		for f = 1:2
			for u = 1:Nunits
				UnitData(u).Net_mean(f) = mean(UnitData(u).NetCounts{f});
				UnitData(u).Net_std(f) = std(UnitData(u).NetCounts{f});
			end
		end

		%--------------------------------------------------------------------------
		%--------------------------------------------------------------------------
		% save data for this set of files
		%--------------------------------------------------------------------------
		%--------------------------------------------------------------------------
		clear tmp
		tmp.UnitIndex = UnitIndex;
		tmp.UnitData = UnitData;
		tmp.BBNinfo = BBN.D.Info;
		tmp.LFHinfo = LFH.D.Info;
		tmp.BBN.ntrials = BBN.ntrials;
		tmp.LFH.ntrials = LFH.ntrials;
		if length(minCommonAtten) == 1
			tmp.minAtten = minCommonAtten * [1 1];
		else
			tmp.minAtten = minCommonAtten;
		end
		tmp.minAttenIndices = minCommonAttenIndices;

		Data{fIndx} = tmp;
		clear tmp;

	end	%%%%%% END OF FILE LOOP

	% write data to mat file
	save(fullfile(outpath, outmatfile), 'Data', '-MAT');
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Perform count of spikes in spike window
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for fIndx = 1:length(Data)
	for u = 1:length(Data{fIndx}.UnitData)
		tmpspikes = Data{fIndx}.UnitData(u).spikes;
		%------------------------------------------------------------------------
		% Count number of spikes in spikeCountWindow{}
		%------------------------------------------------------------------------
		for trial = 1:Data{fIndx}.BBN.ntrials
			for win = 1:Nwin
				spikewin = spikeCountWindow{win};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < tmpspikes{1}{trial}) & ...
										(tmpspikes{1}{trial} < spikewin(2));

				% count them and store in count matrix
				Data{fIndx}.UnitData(u).counts{1}(win, trial) = sum(validSpikes);
			end
		end

		for trial = 1:Data{fIndx}.LFH.ntrials
			for win = 1:Nwin
				spikewin = spikeCountWindow{win};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < tmpspikes{2}{trial}) & ...
										(tmpspikes{2}{trial} < spikewin(2));

				% count them and store in count matrix
				Data{fIndx}.UnitData(u).counts{2}(win, trial) = sum(validSpikes);
			end
		end
		
		for s = 1:2
			for w = 1:length(spikeCountWindow)
				Data{fIndx}.UnitData(u).Count_mean{s}(w) = mean(Data{fIndx}.UnitData(u).counts{s}(w, :));
				Data{fIndx}.UnitData(u).Count_std{s}(w) = std(Data{fIndx}.UnitData(u).counts{s}(w, :));
			end
		end
		
	end
end








%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% export data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% write to .csv file
%--------------------------------------------------------------------------
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
for f = 1:length(Data)
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


