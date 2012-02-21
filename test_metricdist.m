%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load(fullfile(inpath, 'statsdata.mat'));
load('statsdata.mat')

%------------------------------------------------------------------------
% analysis parameters
%------------------------------------------------------------------------
% power transformation for averaging, usually -2
% not used in original SpikeTrain.m function - added back here 
% using jmatrix.m as model
z = -2;


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze LFH Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

Nunits = length(Mlfhstruct{1});
Nconditions = length(Mlfhstruct);
ConditionList = 1:Nconditions


CostVals = [0, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1];

opts.shift_cost = CostVals;
opts.metric_family = 0;
opts.parallel = 0;

uIndex = 1;
% assign data to M struct for simplicity
for c = 1:Nconditions
	M(c) = Mlfhstruct{c}(uIndex);
end

fprintf('File %s, uIndex: %d\n', M(1).file, uIndex);
fprintf('  Probe %d, Cluster %d\n', M(1).probe, M(1).cluster);


Ntrials = zeros(Nconditions, 1);
for c = 1:Nconditions
	Ntrials(c) = M(c).ntrials;
end
TotalTrials = sum(Ntrials);
Times = cell(TotalTrials, 1);
Labels = cell(TotalTrials, 1);
Categories = zeros(TotalTrials, 1);

tindex = 0;
for c = 1:Nconditions
	for n = 1:Ntrials(c)
		tindex = tindex + 1;
		Times{tindex} = M(c).spikes{n};
		Labels{tindex} = ones(size(M(c).spikes{n}));
		Categories(tindex) = c;
	end
end

D = metricdist(1, Times, Labels, opts);

[Dclust, H] = metric_cluster(D, CostVals, Categories, z);


