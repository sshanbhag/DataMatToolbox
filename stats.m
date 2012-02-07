%% statistics

% mnemonics for conditions
PRE = 1;
MILD= 2;
CAT = 3;

% load data if necessary
if ~exist('bbnS', 'var')
	inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
	load(fullfile(inpath, 'BBNrate.mat'), 'bbnS')
end
if ~exist('lfhS', 'var')
	inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
	load(fullfile(inpath, 'LFHrate.mat'), 'lfhS', 'spikeCountWindow', 'bg_spikeCountWindow')
end

% some lengths
Nwin = length(spikeCountWindow);

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



for w = 1:Nbgwin
 %% BBN
 C = [










