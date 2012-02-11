% 1		file
% 2		condition
% 3		unit
% 4		probe
% 5		cluster
% 6		BGspikes
% 7		atten
% 8		ntrials
% 9		spikes
% 10		validspikes
% 11		spikecount
% 12		spikerate
% 13		BG_spiketimes
% 14		BG_spikecount
% 15		BG_spikerate

BGSWEEPTIME_MS = 50;

Windows = {	[0 50] ; ...
				[50 100] ;	...
				[100 150]	;	...
				[150 200]	;	...
			};
WindowSizes = 50;

load statsdata.mat


% build background data
BGcount = cell(3, 1);
BGrate = cell(3, 1);
for c = 1:3
	M = Mbbn{c};
	
	[Nunits(c), Nfields(c)] = size(M);
	tmpcount = cell(Nunits(c), 1);
	tmprate = cell(Nunits(c), 1);
	for u = 1:Nunits(c)
		ntrials = M{u, 8};
		nbgvalues = length(Mbbn{c}{u, 14});
		tmpi = randi(nbgvalues, [ntrials 1]);
		tmpcount{u} = Mbbn{c}{u, 14}(tmpi);
		tmprate{u} = tmpcount{u} ./ (0.001 * BGSWEEPTIME_MS);
	end
	BGcount{c} = tmpcount;
	BGrate{c} = tmprate;
end

% get data for the different windows
Count = cell(length(Windows), 3);
Rate = cell(length(Windows), 3);

for w = 1:length(Windows);

	for c = 1:3
		M = Mbbn{c};

		tmpcount = cell(Nunits(c), 1);
		tmprate = cell(Nunits(c), 1);

		for u = 1:Nunits(c)
			ntrials = M{u, 8};

			tmpcount{u} = M{u, 11}{w};
			tmprate{u} = tmpcount{u} ./ (0.001 * WindowSizes);
		end

		Count{w, c} = tmpcount;
		Rate{w, c} = tmprate;
	end
	
end





%{
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
%}