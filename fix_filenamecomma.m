% fixes issue with comma in filename marker of 
% file '827_12-18-2012--2954_syllables block (30, 40)_Sorted.mat'
close all
clear all
clear classes
%------------------------------------------------------------
%% file names
%------------------------------------------------------------
% location of converted DDF files
inputpath = 'F:\Work\Data\DataWave\batmat\convertedDDF';
matname = '827_12-18-2012--2954_syllables block (30, 40)_Sorted.mat';
% load original mat file
load(fullfile(inputpath, matname));

if exist(fullfile(inputpath, [matname '.ORIG']), 'file')
	disp('File Exists')
	return
end
%------------------------------------------------------------
%% the filename has a comma in it which fucks up 
% reading and processing of marker info.  so, elimnate
% it from the marker and reassign
%------------------------------------------------------------

% loop through event 2 (R channel markers)
for n = 1:D.Event(2).EventCount
	% make copy of event data (marker string, comma delimited)
	tmp = D.Event(2).Data{n};
	% scan CSV string
	tmpR = csvscan(tmp, 0);
	% concatenate the last two elements (separated due to undesired comma)
	ts = [tmpR{end-1} ' ' tmpR{end}];
	% truncate and replace final element by new string
	tmpR = tmpR(1:(end-1));
	tmpR{end} = ts;
	% rebuild comma delimited string
	tmpout = tmpR{1};
	for c = 2:length(tmpR)
		tmpout = [tmpout ',' tmpR{c}];
	end
	% replace marker string in event data
	D.Event(2).Data{n} = tmpout;
end

%------------------------------------------------------------
%% store new struct D
%------------------------------------------------------------

% first, make a copy of the old file
copyfile(fullfile(inputpath, matname), fullfile(inputpath, [matname '.ORIG']));

% then save new D struct
save(fullfile(inputpath, matname), '-MAT', 'D');
