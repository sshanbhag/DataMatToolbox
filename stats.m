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
BGSWEEPTIME_MS = 1000;
Nbgsweeps = round(BGTIME_MS / BGSWEEPTIME_MS);
bgwin = cell(Nbgsweeps, 1);
for b = 1:Nbgsweeps
	bgwin{b} = [b-1 b] * BGSWEEPTIME_MS;
end

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

Nwin = length(spikeCountWindow);
Nbgwin = length(bg_spikeCountWindow);
Nconditions = 3;
ntrials = zeros(Nconditions, 1);
tind = cell(Nconditions, 1);

build_stat_arrays;

return


B = matrixify_data(bbnS);
L = matrixify_data(lfhS);


%% ANOVA, 1 way for each condition


for w = 1:Nwin

	%% BBN 
	C = [B.Count_mean(:, w,  1) B.Count_mean(:, w,  2) B.Count_mean(:, w,  3)];
	clabels = {'Pre', 'Mild', 'Cat'};
	[B.p{w}, B.table{w}, B.stats{w}] = anova1(C, clabels, 'on');
	ylabel('# spikes in window')
	title(sprintf('BBN, Window: %d - %d ms', spikeCountWindow{w}(1), spikeCountWindow{w}(2)))
	figure
	[c,m,h,nms] = multcompare(B.stats{w})
	B.multcomp{w} = struct('c', c, 'm', m, 'h', h, 'names', nms);
	
	%% LFH
	C = [L.Count_mean(:, w,  1) L.Count_mean(:, w,  2) L.Count_mean(:, w,  3)];
	clabels = {'Pre', 'Mild', 'Cat'};
	[L.p{w}, L.table{w}, L.stats{w}] = anova1(C, clabels, 'on');
	ylabel('# spikes in window')
	title(sprintf('LFH, Window: %d - %d ms', spikeCountWindow{w}(1), spikeCountWindow{w}(2)))
	figure
	[c,m,h,nms] = multcompare(L.stats{w})
	L.multcomp{w} = struct('c', c, 'm', m, 'h', h, 'names', nms);

end


%% deal with background









