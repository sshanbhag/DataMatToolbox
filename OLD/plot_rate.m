%% Rate/Condition Plots



if ~exist('bbnS', 'var')
	inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
	load(fullfile(inpath, 'BBNrate.mat'), 'bbnS')
end
if ~exist('lfhS', 'var')
	inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
	load(fullfile(inpath, 'LFHrate.mat'), 'lfhS', 'spikeCountWindow', 'bg_spikeCountWindow')
end

Nwin = length(spikeCountWindow);

B = matrixify_data(bbnS);
L = matrixify_data(lfhS);

PRE = 1;
MILD= 2;
CAT = 3;

% Scatter plots


for w = 1:Nwin

	%% BBN 

	figure(w)

	% plot PRE vs MILD
	subplot(2,4,1)
	plot(B.Count_mean(:, w,  PRE), B.Count_mean(:, w,  MILD), '.')
	mx = max([B.Count_mean(:, w,  PRE); B.Count_mean(:, w,  MILD)]);
	xlim([0 mx]);
	ylim([0 mx]);
	axis square
	line([0 mx], [0 mx], 'Color', [1 0 1])
	xlabel('Pre (1)')
	ylabel('Mild (2)')
	title(sprintf('BBN, Window = [%d %d] ms', spikeCountWindow{w}(1), spikeCountWindow{w}(2)))

	% plot PRE vs. CAT
	subplot(2,4,2)
	plot(B.Count_mean(:, w,  PRE), B.Count_mean(:, w,  CAT), '.')
	mx = max([B.Count_mean(:, w,  PRE); B.Count_mean(:, w,  CAT)]);
	xlim([0 mx]);
	ylim([0 mx]);
	axis square
	line([0 mx], [0 mx], 'Color', [1 0 1])
	xlabel('Pre (1)')
	ylabel('Cat (3)')

	% plot MILD vs. CAT
	subplot(2,4,3)
	plot(B.Count_mean(:, w,  MILD), B.Count_mean(:, w,  CAT), '.')
	mx = max([B.Count_mean(:, w,  MILD); B.Count_mean(:, w,  CAT)]);
	xlim([0 mx]);
	ylim([0 mx]);
	axis square
	line([0 mx], [0 mx], 'Color', [1 0 1])
	xlabel('Mild (2)')
	ylabel('Cat (3)')

	subplot(2,4,4)
	C = [B.Count_mean(:, w,  1) B.Count_mean(:, w,  2) B.Count_mean(:, w,  3)];
	clabels = {'Pre', 'Mild', 'Cat'};
	boxplot(C, clabels)
	ylabel('# spikes in window')
	title(sprintf('BBN, Window: %d - %d ms', spikeCountWindow{w}(1), spikeCountWindow{w}(2)))

	%% LFH
	% plot PRE vs MILD
	subplot(2,4,5)
	plot(L.Count_mean(:, w,  PRE), L.Count_mean(:, w,  MILD), '.')
	mx = max([L.Count_mean(:, w,  PRE); L.Count_mean(:, w,  MILD)]);
	xlim([0 mx]);
	ylim([0 mx]);
	axis square
	line([0 mx], [0 mx], 'Color', [1 0 1])
	xlabel('Pre (1)')
	ylabel('Mild (2)')
	title(sprintf('LFH, Window = [%d %d] ms', spikeCountWindow{w}(1), spikeCountWindow{w}(2)))

	% plot PRE vs. CAT
	subplot(2,4,6)
	plot(L.Count_mean(:, w,  PRE), L.Count_mean(:, w,  CAT), '.')
	mx = max([L.Count_mean(:, w,  PRE); L.Count_mean(:, w,  CAT)]);
	xlim([0 mx]);
	ylim([0 mx]);
	axis square
	line([0 mx], [0 mx], 'Color', [1 0 1])
	xlabel('Pre (1)')
	ylabel('Cat (3)')

	% plot MILD vs. CAT
	subplot(2,4,7)
	plot(L.Count_mean(:, w,  MILD), L.Count_mean(:, w,  CAT), '.')
	mx = max([L.Count_mean(:, w,  MILD); L.Count_mean(:, w,  CAT)]);
	xlim([0 mx]);
	ylim([0 mx]);
	axis square
	line([0 mx], [0 mx], 'Color', [1 0 1])
	xlabel('Mild (2)')
	ylabel('Cat (3)')

	subplot(2,4,8)
	C = [L.Count_mean(:, w,  PRE) L.Count_mean(:, w,  MILD) L.Count_mean(:, w,  CAT)];
	clabels = {'Pre', 'Mild', 'Cat'};
	boxplot(C, clabels)
	ylabel('# spikes in window')
	title(sprintf('LFH, Window = [%d %d] ms', spikeCountWindow{w}(1), spikeCountWindow{w}(2)))
	ylim([0 mx])
end















