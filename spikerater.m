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
% define some options
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles';

types_to_process = {'BBN', 'LFH'};

outpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles';

% size of time window after stimulus onset in which to count spikes
RESPONSEWINDOW_MS = 100;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some constants/settings
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% default sound onset time time
soundOnsetTime = 0;
spikeCountWindow{1} = [0 RESPONSEWINDOW_MS];
Nwin = length(spikeCountWindow);

%--------------------------------------------------------------------------
% define some plot constants/settings
%--------------------------------------------------------------------------

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
%--------------------------------------------------------------------------
Data = cell(1, match_index);

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
	if BBN.D.Info.Nunits ~= LFH.D.Info.Nunits
		warning('%s: # Units mismatch!!!!!!!!!!!\n', mfilename)
		fprintf('\tBBN file has %d\tunits\n', BBN.D.Info.Nunits);
		fprintf('\tLFH file has %d\tunits\n\n', LFH.D.Info.Nunits);
		
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
		% and store unit numbers in UnitIndex[BBN.D.Info.UnitList(:, 3) LFH.D.Info.UnitList(:, 3)]
		if BBN.D.Info.Nunits > LFH.D.Info.Nunits
			UnitList = zeros(LFH.D.Info.Nunits, 2);
			UnitIndex = UnitList;
			u = 0;
			for p = 1:BBN.D.Info.Nunits
				lFlag = 0;
				for l = 1:LFH.D.Info.Nunits
					if all(blist(p, 1:2) == llist(l, 1:2))
						lFlag = l;
					end
				end
				
				if lFlag
					u = u + 1;
					UnitList(u, :) = blist(p, 1:2);
					UnitIndex(u, :) = [p lFlag];
				end
			end
		else
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
	nVarA = length(BBN.Stimulus.Var);
	attenindex = 0;
	for n = 1:nVarA
		if strcmp(BBN.Stimulus.Var(n).name, 'AttenuationR')
			attenindex = n;
		end
	end
	if attenindex
		attenValsA = BBN.Stimulus.Var(attenindex).values;
	else
		error('%s: attenuation not found in Var list for file %s', mfilename, bbnfile);
	end
	NattenA = length(attenValsA);
	
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
	for n = 1:NattenA
		cind = find(attenValsA(n) == attenValsB);
		if ~isempty(cind)
			ccount = ccount + 1;
			commonAtten(ccount) = attenValsA(n);
		end
	end
	if ccount
		minCommonAtten = min(commonAtten);
		minCommonAttenIndices = ...
			[find(attenValsA == minCommonAtten) find(attenValsB == minCommonAtten)];
	else
		warning('%s: no common attenuation values in files %s and %s', mfilename, bbnfile, lfhfile);
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
		% Count number of spikes in spikeCountWindow{}
		%------------------------------------------------------------------------
		for trial = 1:BBN.ntrials
			for win = 1:Nwin
				spikewin = spikeCountWindow{win};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < UnitData(u).spikes{1}{trial}) & ...
										(UnitData(u).spikes{1}{trial} < spikewin(2));

				% count them and store in count matrix
				UnitData(u).counts{1}(win, trial) = sum(validSpikes);
			end
		end
		
		for trial = 1:LFH.ntrials
			for win = 1:Nwin
				spikewin = spikeCountWindow{win};
				% find spikes that are within the current window
				validSpikes = (spikewin(1) < UnitData(u).spikes{2}{trial}) & ...
										(UnitData(u).spikes{2}{trial} < spikewin(2));

				% count them and store in count matrix
				UnitData(u).counts{2}(win, trial) = sum(validSpikes);
			end
		end

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

			for w = 1:length(spikeCountWindow)
				UnitData(u).Count_mean{f}(w) = mean(UnitData(u).counts{f}(w, :));
				UnitData(u).Count_std{f}(w) = std(UnitData(u).counts{f}(w, :));
			end

		end
	end
	
	
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	% save data for this set of files
	%--------------------------------------------------------------------------
	%--------------------------------------------------------------------------
	clear tmp
	tmp.UnitData = UnitData;
	tmp.BBNinfo = BBN.D.Info;
	tmp.LFHinfo = LFH.D.Info;
	Data{fIndx} = tmp;
	clear tmp;
	
end	%%%%%% END OF FILE LOOP

return




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% export data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% write to rate file
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fp = fopen('tmpout.csv', 'w');

fprintf(fp, 'SpikeCountWindowLength,%f,msec,\n', RESPONSEWINDOW_MS);

fprintf(fp, 'File,');
fprintf(fp, 'UnitNum,');
fprintf(fp, 'Probe,');
fprintf(fp, 'Cluster,');
fprintf(fp, 'BGcount,');
for a = 1:Natten
	fprintf(fp, '%d dB Net_mean,', indata.attenVals(a));
end
for a = 1:Natten
	fprintf(fp, '%d dB Net_sd,', indata.attenVals(a));
end
for a = 1:Natten
	fprintf(fp, '%d dB Count_mean,', indata.attenVals(a));
end
for a = 1:Natten
	fprintf(fp, '%d dB Count_sd,', indata.attenVals(a));
end

fprintf(fp, '\n');

for u = 1:Nunits
	fprintf(fp, '%s,', indata.filename{1});
	fprintf(fp, '%d,', UnitData(u).UnitInfo.unit);
	fprintf(fp, '%d,', UnitData(u).UnitInfo.probe);
	fprintf(fp, '%d,', UnitData(u).UnitInfo.cluster);
	fprintf(fp, '%f,', UnitData(u).BG_count);
	
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Net_mean(a));
	end
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Net_std(a));
	end
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Count_mean{1, a});
	end
	for a = 1:Natten
		fprintf(fp, '%f,', UnitData(u).Count_std{1, a});
	end
	
	
	
	fprintf(fp, '\n');
end

if fp ~= 1
	fclose(fp);
end


