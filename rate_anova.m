%-------------------------------------------------------------------
% rate_anova
%-------------------------------------------------------------------



%% -----------------------------------------------------------------
%-------------------------------------------------------------------
% Some Constants
%-------------------------------------------------------------------

%-------------------------------------------------------------------
%-------------------------------------------------------------------
% mnemonics for conditions
%-------------------------------------------------------------------
PRE = 1;
MILD= 2;
CAT = 3;

load statsdata.mat

Nconditions = length(Mbbnstruct);
Nunits = length(Mbbnstruct{1});

for u = 1:1
	
	% make a tmp copy of current unit's data
	for c = 1:3
		M(c) = Mbbnstruct{c}(u);
	end
	
	C = 
	




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