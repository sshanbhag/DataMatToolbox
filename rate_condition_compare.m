%------------------------------------------------------------------------
% rate_condition_compare.m
%------------------------------------------------------------------------
% script for spike rate computation and analysis from DataMat toolbox
%------------------------------------------------------------------------
% Some Notes:
%------------------------------------------------------------------------

% clear all

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some constants/settings
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% size of time window after stimulus onset in which to count spikes
RESPONSEWINDOW_MS = 100;

% default sound onset time time
soundOnsetTime = 0;
spikeCountWindow{1} = [0 100];
spikeCountWindow{2} = [0 200];
spikeCountWindow{3} = [0 300];

Nwin = length(spikeCountWindow);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some options
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles';

types_to_process = {'BBN', 'LFH'};

outpath = [inpath '/output'];
matrootname = 'Data_BBN_LFH'
outmatfile = [matrootname '.mat'];
outcountmatfile = sprintf('%s_%dms.mat', matrootname, RESPONSEWINDOW_MS);
outcsvfile = sprintf('%s_%dms.csv', matrootname, RESPONSEWINDOW_MS);

if ~exist(outpath, 'dir')
	mkdir(outpath);
end

%------------------------------------------------------------------------
% List of files and units to analyze
%------------------------------------------------------------------------
% 1) need to find matching pairs of files
% get all BBN files first, then search for matching LFH files?
%------------------------------------------------------------------------

% get list of files within the inpath directory 
% and store info in fstruct array
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

%{
tmpfiles = cell(bbnCount, 1);
for n = 1:bbnCount
	tmpfiles{n} = fullfile(inpath, bbnFiles{n});
end
bbnData = spikerater(tmpfiles, spikeCountWindow);

tmpfiles = cell(lfhCount, 1);
for n = 1:lfhCount
	tmpfiles{n} = fullfile(inpath, lfhFiles{n});
end
lfhData = spikerater(tmpfiles, spikeCountWindow);

save(fullfile(outpath, 'data.mat'), 'bbnData', 'lfhData', '-MAT')
%}

if ~exist('bbnData', 'var')
	fprintf('Loading data file %s\n', fullfile(outpath, 'data.mat'));
	load(fullfile(outpath, 'data.mat'), 'bbnData', 'lfhData', '-MAT')
end

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
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% to search for files for same experiment but different conditions, 
% need to remove _xxx characters and condition id
% store these is bbn_idstr and lfh_idstr
%--------------------------------------------------------------------------

% go through the list of BBN filenames and remove the _BBN chars and 
% last 2 characters that correspond to '_<condition>'
for n = 1:bbnCount
	% parse file name
	[~, tmpfname, ext] = fileparts(bbnFiles{n});
	% remove _BBN
	tmpfname = regexprep(tmpfname, ('_BBN'), '');
	% remove last 2 chars and store in bbn_idstr cell array
	bbn_idstr{n} = tmpfname(1:(end - 2));
end
% Repeat this for LFH files
for n = 1:lfhCount
	% parse file name
	[~, tmpfname, ext] = fileparts(lfhFiles{n});
	% remove _LFH
	tmpfname = regexprep(tmpfname, ('_LFH'), '');
	% remove last 2 chars and store in lfh_idstr cell array
	lfh_idstr{n} = tmpfname(1:(end - 2));
end
	
%--------------------------------------------------------------------------
% find matching filenames/conditions
%--------------------------------------------------------------------------
bbnList = find_condition_instances(bbn_idstr, bbnFiles);
lfhList = find_condition_instances(lfh_idstr, lfhFiles);

% get the conditions and store in bbnList
for n = 1:length(bbnList)
	tmp = zeros(1, length(bbnList{n, 1}));
	for c = 1:length(bbnList{n, 1})
	 tmp(c) = bbnData{bbnList{n, 1}(c)}.Info.condition;
	end
	bbnList{n, 4} = tmp;
end
for n = 1:length(lfhList)
	tmp = zeros(1, length(lfhList{n, 1}));
	for c = 1:length(lfhList{n, 1})
	 tmp(c) = lfhData{lfhList{n, 1}(c)}.Info.condition;
	end
	lfhList{n, 4} = tmp;
end

% pre-allocate Data cell array to hold data for each pair of BBN and LFH files
% BBNData = cell(1, length(bbnList));

%--------------------------------------------------------------------------
% limit list to analyze to those files for which there are all three
% conditions (1, 2, 3)
%--------------------------------------------------------------------------
% valid___List{}
% 	{n, 1}	->	index into bbnData and bbnFiles
% 	{n, 1}	->	file search strings (no _BBN or _LFH or condition value)
% 	{n, 2}	->	full .mat file name
% 	{n, 4}	->	condition #
%--------------------------------------------------------------------------

% BBN files
n = 0;
for bIndx = 1:length(bbnList)
	clist = bbnList{bIndx, 4};
	if (~isempty(clist)) && (length(clist) >= 3)
		% see if desired conditions are in the list
		ctest = zeros(3, length(clist));
		ctest(1, :) = (clist == 1);
		ctest(2, :) = (clist == 2);
		ctest(3, :) = (clist == 3);
		if sum(sum(ctest)) == 3
			for c = 1:3
				cind(c) = find(ctest(c, :));
			end

			% if so, sort and store the file information
			n = n + 1;
			validBBNList(n, :) = bbnList(bIndx, :);
			for c = 1:3
				for j = 1:4
					validBBNList{n, j} = validBBNList{n, j}(cind);
				end
			end
		end
	end
end

% LFH files
n = 0;
for lIndx = 1:length(lfhList)
	clist = lfhList{lIndx, 4};
	if (~isempty(clist)) && (length(clist) >= 3)
		% see if desired conditions are in the list
		ctest = zeros(3, length(clist));
		ctest(1, :) = (clist == 1);
		ctest(2, :) = (clist == 2);
		ctest(3, :) = (clist == 3);
		if sum(sum(ctest)) == 3
			for c = 1:3
				cind(c) = find(ctest(c, :));
			end

			% if so, sort and store the file information
			n = n + 1;
			validLFHList(n, :) = lfhList(lIndx, :);
			for c = 1:3
				for j = 1:4
					validLFHList{n, j} = validLFHList{n, j}(cind);
				end
			end
		end
	end
end



% loop through file list
for fIndx = 1:length(validBBNList)
	bfiles = validBBNList{fIndx, 3};
	bindex = validBBNList{fIndx, 1};

	if isempty(bfiles)
		error('%s: empty bfiles', mfilename)
	end
	
	Ncond = length(bindex);

	% load BBN file data into list
	BBN = cell(1, Ncond);
	Blist = cell(1, Ncond);
	for n = 1:Ncond
		BBN{n} = bbnData{bindex(n)};
		% columns of blist arrays are probe #, cluster # and unit #
		Blist{n} = BBN{n}.Info.UnitList;
	end
	
	[unitM, matches] = match_units(Blist);
	if isempty(unitM)
		warning('%s: no unit matches found', mfilename)
	end
	validBBNList{fIndx, 5} = matches;
	validBBNList{fIndx, 6} = unitM;
end
	

	%--------------------------------------------------------------------------
	% determine number of Attenuation levels by looking at the Stimulus.Var
	% information in each condition's data structure
	%--------------------------------------------------------------------------
for fIndx = 1:length(validBBNList)
	bfiles = validBBNList{fIndx, 3};
	bindex = validBBNList{fIndx, 1};

	if isempty(bfiles)
		error('%s: empty bfiles', mfilename)
	end
	
	Ncond = length(bindex);

	% load BBN file data into list
	Blist = cell(1, Ncond);
	for c = 1:Ncond
		BBN{c} = bbnData{bindex(c)};
		% columns of Blist arrays are probe #, cluster # and unit #
		Blist{c} = BBN{c}.Stimulus.Var;
		% # of attenuation values in BBN file
		nVar = length(Blist{c});

		attenindex = 0;
		for n = 1:nVar
			% find R channel attenuation
			if strcmp(Blist{c}(n).name, 'AttenuationR')
				attenindex = n;
			end
		end
		
		% if found, store the values
		if attenindex
			attenVals{c} = Blist{c}(attenindex).values;
		else
			% if not, error
			error('%s: attenuation not found in Var list for file %s', mfilename, bfiles{c});
		end
		Natten(c) = length(attenVals{c});				
	end

	[m, l] = match_atten(attenVals);
	validBBNList{fIndx, 7} = l;
	validBBNList{fIndx, 8} = m;
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% export data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% write to .csv file
%--------------------------------------------------------------------------

fp = 1;

% fp = fopen(fullfile(outpath, 'temp.csv'), 'w');

fprintf(fp, 'SpikeCountWindowLength,%f,msec,\n', RESPONSEWINDOW_MS);

fprintf(fp, 'File,');
fprintf(fp, 'UnitNum,');
fprintf(fp, 'Probe,');
fprintf(fp, 'Cluster,');
fprintf(fp, 'BGcount,');
fprintf(fp, 'Atten_dB,');
fprintf(fp, 'Net_mean,');
fprintf(fp, 'Net_sd,');
fprintf(fp, 'Count_mean,');
fprintf(fp, 'Count_sd,');
fprintf(fp, '\n');

rindx = 0;
for f = 1:length(validBBNList)
	
	validDataIndices = validBBNList{f, 1};
	validUnitIndices = validBBNList{f, 5};
	validAttenIndices = validBBNList{f, 7};
	
	Nconditions = length(validDataIndices)
	[Nunits, tmp] = size(validUnitIndices);

	for c = 1:Nconditions
		tmpD = bbnData{validDataIndices(c)}
		attIndex = validAttenIndices(1, c);

		for u = 1:Nunits
			% columns in validUnitIndices correspond to conditions/files
			vUI = validUnitIndices(u, c);
			tmpU = tmpD.UnitData(vUI);
keyboard	
			fprintf(fp, '%s,', tmpD.Info.file);
			fprintf(fp, '%d,', tmpU.UnitInfo.unit);
			fprintf(fp, '%d,', tmpU.UnitInfo.probe);
			fprintf(fp, '%d,', tmpU.UnitInfo.cluster);
			fprintf(fp, '%.4f,', tmpU.BG_count);
			fprintf(fp, '%.4f,', tmpD.AttenVals(attIndex));
			fprintf(fp, '%.4f,', tmpU.Net_mean(attIndex));
			fprintf(fp, '%.4f,', tmpU.Net_std(attIndex));
			
			% loop through spike count windows
			for w = 1:Nwin
				fprintf(fp, '%.4f,', tmpU.Count_mean{attIndex}(w));
				fprintf(fp, '%.4f,', tmpU.Count_std{attIndex}(w));
			end
			fprintf(fp, '\n');
		
		rindx = rindx + 1;
		col = 1;
		Rdata{rindx, col} = tmpD.BBNinfo.file;	col = col + 1;
		Rdata{rindx, col} = tmpU.UnitInfo.unit;	col = col + 1;
		Rdata{rindx, col} = tmpU.UnitInfo.probe;	col = col + 1;
		Rdata{rindx, col} = tmpU.UnitInfo.cluster;	col = col + 1;
		Rdata{rindx, col} = tmpU.BG_count;	col = col + 1;
		Rdata{rindx, col} = tmpD.AttenVals(attIndex);	col = col + 1;
		Rdata{rindx, col} = tmpU.Net_mean(attIndex);	col = col + 1;
		Rdata{rindx, col} = tmpU.Net_std(attIndex);	col = col + 1;
		Rdata{rindx, col} = tmpU.Count_mean{1}(attIndex);	col = col + 1;
		Rdata{rindx, col} = tmpU.Count_std{1}(attIndex);	col = col + 1;
	end
end
end
if fp ~= 1
	fclose(fp);
end

