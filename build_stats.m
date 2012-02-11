%-------------------------------------------------------------------
% stats
%-------------------------------------------------------------------

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
% mnemonics for conditions
%-------------------------------------------------------------------
PRE = 1;
MILD= 2;
CAT = 3;
%-------------------------------------------------------------------
% colors for different conditions
%-------------------------------------------------------------------
tcolors{PRE}	= 'b';
tcolors{MILD}	= 'g';
tcolors{CAT}	= 'r';
%-------------------------------------------------------------------
% background data parameters
%-------------------------------------------------------------------
BGTIME_MS = 5000;
BGSWEEPTIME_MS = 50;

%-------------------------------------------------------------------
% file paths
%-------------------------------------------------------------------
inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
outpath = [inpath '/stats'];
if ~exist(outpath, 'dir')
	mkdir(outpath);
end


%% -----------------------------------------------------------------
%-------------------------------------------------------------------
% Load data if necessary
%-------------------------------------------------------------------
%-------------------------------------------------------------------
if ~exist('bbnData', 'var')
	load(fullfile(inpath, 'BBNrate.mat'), 'validBBNList', 'bbnData', ...
						'bbnS', 'spikeCountWindow', 'bg_spikeCountWindow')
end
if ~exist('lfhData', 'var')
	load(fullfile(inpath, 'LFHrate.mat'),  'validLFHList', 'lfhData', ...
						'lfhS', 'spikeCountWindow', 'bg_spikeCountWindow')
end

Windows = {	[0 50] ; ...
				[50 100] ;	...
				[100 150]	;	...
				[150 200]	;	...
			};
BGWindow = 50;
Clist = [1 2 3];

[Mbbn, Mbbnstruct, Mfields] = build_stat_arrays(bbnData, validBBNList, Clist, Windows, BGSWEEPTIME_MS, BGTIME_MS);
[Mlfh, Mlfhstruct] = build_stat_arrays(lfhData, validLFHList, Clist, Windows, BGSWEEPTIME_MS, BGTIME_MS);

save(	'statsdata.mat', 'Mbbn', 'Mbbnstruct', 'Mfields', ...
		'Mlfh', 'Mlfhstruct', 'Windows', 'BGWindow', 'Clist', '-MAT');





