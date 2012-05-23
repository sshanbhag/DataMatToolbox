
if ~exist('Dbbn', 'var')
	load spike_distance_output.mat
end


% plot, for each unit, the information -vs- cost
for u = 1:length(Dbbn)
	figure(1)
	hval = zeros(length(Dbbn(u).H), 1);
	for n = 1:length(Dbbn(u).H)
		hval(n) = Dbbn(u).H(n).value;
	end
	plot(opts.shift_cost, hval, '.-')
	xlabel('cost/msec')
	ylabel('H (bits)')
	[a, aind] = max(hval);
	hold on
		plot(CostVals(aind), a, 'ro')
	hold off
	title({	sprintf('Dataset: %s', Mbbnstruct{1}(u).file),	...
				sprintf('Probe: %d Cluster: %d', Mbbnstruct{1}(u).probe, Mbbnstruct{1}(u).cluster),	...
				sprintf('qmax=%.2f', CostVals(aind))	}, ...
				'Interpreter', 'none');
	MaxH_BBN(u) = a;
	MaxC_BBN(u) = opts.shift_cost(aind);
end

for u = 1:length(Dlfh)
	figure(1)
	hval = zeros(length(Dlfh(u).H), 1);
	for n = 1:length(Dlfh(u).H)
		hval(n) = Dlfh(u).H(n).value;
	end
	plot(opts.shift_cost, hval, '.-')
	xlabel('cost/msec')
	ylabel('H (bits)')
	[a, aind] = max(hval);
	hold on
		plot(CostVals(aind), a, 'ro')
	hold off
	title({	sprintf('Dataset: %s', Mlfhstruct{1}(u).file),	...
				sprintf('Probe: %d Cluster: %d', Mlfhstruct{1}(u).probe, Mlfhstruct{1}(u).cluster),	...
				sprintf('qmax=%.2f', CostVals(aind))	}, ...
				'Interpreter', 'none');
	MaxH_LFH(u) = a;
	MaxC_LFH(u) = opts.shift_cost(aind);
end



figure
maxvalH = max([MaxH_BBN MaxH_LFH]);
subplot(221)
boxplot(MaxH_BBN);
ylim([0 maxvalH]);
subplot(222)
boxplot(MaxH_LFH)
ylim([0 maxvalH]);

maxvalC = max([MaxC_BBN MaxC_LFH]);
subplot(223)
boxplot(MaxC_BBN);
ylim([0 maxvalC]);
subplot(224)
boxplot(MaxC_LFH)
ylim([0 maxvalC]);

figure
subplot(121)
plot(MaxC_BBN, MaxH_BBN ,'.');
xlim([0 maxvalC])
ylim([0 maxvalH])

subplot(122)
plot(MaxC_LFH, MaxH_LFH, '.');
xlim([0 maxvalC])
ylim([0 maxvalH])





