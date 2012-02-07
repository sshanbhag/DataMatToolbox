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

% spike count windows
spikeCountWindow{1} = [0 100];
spikeCountWindow{2} = [0 200];
spikeCountWindow{3} = [0 300];
spikeCountWindow{4} = [0 50];
spikeCountWindow{5} = [50 100];
spikeCountWindow{6} = [100 300];
Nwin = length(spikeCountWindow);

% background count windows
bg_spikeCountWindow(1) = 50;
bg_spikeCountWindow(2) = 100;
bg_spikeCountWindow(3) = 200;
bg_spikeCountWindow(4) = 300;
Nbgwin = length(bg_spikeCountWindow);

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

	% check if this is one of the types to process - 
	% if so, add to appropriate list
	if strfind(fname, 'BBN')
		bbnCount = bbnCount + 1;
		bbnFiles{bbnCount} = fname;
	elseif strfind(fname, 'LFH')
		lfhCount = lfhCount + 1;
		lfhFiles{lfhCount} = fname;
	end
end

%% -------------------------------------------------------------------------
% Now compute spike rate for all files and units
%--------------------------------------------------------------------------

tmpfiles = cell(bbnCount, 1);
for n = 1:bbnCount
	tmpfiles{n} = fullfile(inpath, bbnFiles{n});
end
bbnData = spikerater(tmpfiles, spikeCountWindow, bg_spikeCountWindow);

tmpfiles = cell(lfhCount, 1);
for n = 1:lfhCount
	tmpfiles{n} = fullfile(inpath, lfhFiles{n});
end
lfhData = spikerater(tmpfiles, spikeCountWindow, bg_spikeCountWindow);

save(fullfile(outpath, 'data.mat'), 'bbnData', 'lfhData', '-MAT')

% if ~exist('bbnData', 'var')
% 	fprintf('Loading data file %s\n', fullfile(outpath, 'data.mat'));
% 	load(fullfile(outpath, 'data.mat'), 'bbnData', 'lfhData', '-MAT')
% end

%% -------------------------------------------------------------------------
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


%% ------------------------------------------------------------------------
% limit list to analyze to those files for which there are all three
% conditions (1, 2, 3)
%--------------------------------------------------------------------------
% valid___List{}
% 	{n, 1}	->	index into bbnData and bbnFiles
% 	{n, 1}	->	file search strings (no _BBN or _LFH or condition value)
% 	{n, 2}	->	full .mat file name
% 	{n, 4}	->	condition #
%--------------------------------------------------------------------------

validBBNList = find_valid_data(bbnData, bbnList);

validLFHList = find_valid_data(lfhData, lfhList);

%% ------------------------------------------------------------------------
%--------------------------------------------------------------------------
% export data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% write to .csv file
%--------------------------------------------------------------------------
excludeList = [ 0 ];

bbnR = export_data(fullfile(outpath, 'BBN.csv'), ...
							bbnData, ...
							validBBNList, ...
							spikeCountWindow	);

bbnS = export_data_by_condition(fullfile(outpath, 'BBN_bycondition.csv'), ...
							bbnData, ...
							validBBNList, ...
							spikeCountWindow, ...
							excludeList	);
						
save(fullfile(outpath, 'BBNrate.mat'), 'bbnData', 'bbnList', 'validBBNList', ...
									'bbnR', 'bbnS', ...
									'spikeCountWindow', 'bg_spikeCountWindow', '-MAT');
						
lfhR = export_data(fullfile(outpath, 'LFH.csv'), ...
							lfhData, ...
							validLFHList, ...
							spikeCountWindow	);

lfhS = export_data_by_condition(fullfile(outpath, 'LFH_bycondition.csv'), ...
							lfhData, ...
							validLFHList, ...
							spikeCountWindow, ...
							excludeList	);
						
save(fullfile(outpath, 'LFHrate.mat'), 'lfhData', 'lfhList','validLFHList', ...
									'lfhR', 'lfhS', ...
									'spikeCountWindow', 'bg_spikeCountWindow', '-MAT');

