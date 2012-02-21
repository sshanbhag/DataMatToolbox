%------------------------------------------------------------------------
% spike_distance.m
%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%
%	uses code from expfun project by Chris Sumner as well
% 	as bits from the spike train toolbox
%------------------------------------------------------------------------
% Created: 1 December 2011 (SJS)
%
% Revisions:
%	9 Dec 2011 (SJS): vectorized Dave computation, added comments
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear old bits
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear all
close all

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Settings for analysis
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% analysis parameters
%------------------------------------------------------------------------
% power transformation for averaging, usually -2
% not used in original SpikeTrain.m function - added back here 
% using jmatrix.m as model
z = -2;

% Cost of changing time, A cost of 0 becomes rate based, a cost of 2/cost
% makes it cheaper to add and delete spikes than to move them.
% change this value to an appropriate valuse for the data, the hard part is picking
% an appropriate value for all of the data
CostVals = [0, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1];
Ncost = length(CostVals);

% options for metricdist.mex function
opts.shift_cost = CostVals;
opts.metric_family = 0;
opts.parallel = 0;

%------------------------------------------------------------------------
% File Settings
%------------------------------------------------------------------------

inpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles/output';
outpath = [inpath '/stats'];
if ~exist(outpath, 'dir')
	mkdir(outpath);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load(fullfile(inpath, 'statsdata.mat'));
load('statsdata.mat')

%{
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze BBN Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

Nunits = length(Mbbnstruct{1});
Nconditions = length(Mbbnstruct);
ConditionList = 1:Nconditions

%------------------------------------------------------------------------
% loop units
%------------------------------------------------------------------------
for uIndex = 1:Nunits
	
	% assign data to M struct for simplicity
	for condition = 1:Nconditions
		M(condition) = Mbbnstruct{condition}(uIndex);
	end
	fprintf('File %s\n', M(1).file);
	fprintf('  Probe %d, Cluster %d\n', M(1).probe, M(1).cluster);
	
	Ntrials = M(1).ntrials;

	% loop through cost values
	for Cindex = 1:Ncost
		% get current cost value from list
		Cost = CostVals(Cindex);
		fprintf('\tCost = %f\n', Cost);

		% initialize Ddistance array
		Ddistance = zeros(Nconditions, Nconditions);

		% loop for the number of conditions
		for C1 = 1:Nconditions

			% loop for the number of trials;
			for T1 = 1:Ntrials

				% initialize Dspike
				Dspike = zeros(Nconditions, Ntrials);

				% for each pair of conditions (stimuli), compute spike
				% distance metric for all trials (except for this
				% trial/condition pair with itself)
				for C2 = 1:Nconditions
					for T2=1:Ntrials
						Dspike(C2,T2) = (spkd(M(C1).spikes{T1}, M(C2).spikes{T2}, Cost))^z;
					end;
				end;

				% find values of Dspike that are 'Inf' and set them to zero
				[a,b] = find(Dspike==Inf);
				Dspike(a,b) = 0;

				% compute overall sum across trials (rows of Dspike matrix)
				Dsum = sum(Dspike,2);

				% loop through conditions, using appropriate # of trials to 
				% compute average for situation in which loop condition == current 
				% overall condition for this cost value
				Dave = zeros(Nconditions, 1);
				currentCondition = ConditionList(C1 == ConditionList);
				otherNconditions = ConditionList(C1 ~= ConditionList);
				Dave(currentCondition) = (Dsum(currentCondition)./(Ntrials-1)).^(1/z);
				Dave(otherNconditions) = (Dsum(otherNconditions)./Ntrials).^(1/z);

				% find indices of minimum avg distance
				[c ,d] = min(Dave);
				% increment count for this overall condition
				Ddistance(C1, d) = Ddistance(C1, d) + 1;

				% clear some things
				clear Dspike Dsum Dave;
			end
		end

		% compute some stats from Ddistance (confusion) matrix
		[Hvalue, Hrand, Hprob] = Hstat(Ddistance);

		% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
		% xlswrite(ExsName, Ddistance)
		fprintf('\tHvalue: %f\n', Hvalue);
		fprintf('\tHrand: %f\n', Hrand);
		fprintf('\tHprob: %f\n', Hprob);
		fprintf('\tDdistance = \n');
		disp(Ddistance);

		% store data
		Dbbn(uIndex).Ddistance{Cindex} = Ddistance;
		Dbbn(uIndex).Hvalue(Cindex) = Hvalue;
		Dbbn(uIndex).Hrand(Cindex) = Hrand;
		Dbbn(uIndex).Hprob(Cindex) = Hprob;

		clear Ddistance SpkTrains;
	end
end
%}

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze BBN Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

Nunits = length(Mbbnstruct{1});
Nconditions = length(Mbbnstruct);
ConditionList = 1:Nconditions

%------------------------------------------------------------------------
% loop units
%------------------------------------------------------------------------
for uIndex = 1:Nunits
	
	% assign data to M struct for simplicity
	for c = 1:Nconditions
		M(c) = Mbbnstruct{c}(uIndex);
	end
	
	fprintf('File %s, uIndex: %d\n', M(1).file, uIndex);
	fprintf('  Probe %d, Cluster %d\n', M(1).probe, M(1).cluster);
	
	% get the # of trials for each condition
	Ntrials = zeros(Nconditions, 1);
	for c = 1:Nconditions
		Ntrials(c) = M(c).ntrials;
	end
	% TotalTrials is total # of trials across conditions - this
	% will be used to allocate vectors for input into metricdist.mex function
	TotalTrials = sum(Ntrials);
	Times = cell(TotalTrials, 1);
	Labels = cell(TotalTrials, 1);
	Categories = zeros(TotalTrials, 1);

	tindex = 0;
	% loop through conditons
	for c = 1:Nconditions
		% loop through trials for this condition
		for n = 1:Ntrials(c)
			% assign spiketimes and labels (unit #) and category (condition)
			% to vectors used for metricdist.mex function
			tindex = tindex + 1;
			Times{tindex} = M(c).spikes{n};
			Labels{tindex} = uIndex*ones(size(M(c).spikes{n}));
			Categories(tindex) = c;
		end
	end

	% compute spike distance metric
	D = metricdist(1, Times, Labels, opts);

	% compute some stats from Ddistance (confusion) matrix
	[Dclust, H] = metric_cluster(D, CostVals, Categories, z);

	% store data
	Dbbn(uIndex).D = D;
	Dbbn(uIndex).Dclust = Dclust;
	Dbbn(uIndex).H = H;
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze LFH Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

Nunits = length(Mlfhstruct{1});
Nconditions = length(Mlfhstruct);
ConditionList = 1:Nconditions

%------------------------------------------------------------------------
% loop units
%------------------------------------------------------------------------
for uIndex = 1:Nunits
	
	% assign data to M struct for simplicity
	for c = 1:Nconditions
		M(c) = Mlfhstruct{c}(uIndex);
	end
	
	fprintf('File %s, uIndex: %d\n', M(1).file, uIndex);
	fprintf('  Probe %d, Cluster %d\n', M(1).probe, M(1).cluster);
	
	% get the # of trials for each condition
	Ntrials = zeros(Nconditions, 1);
	for c = 1:Nconditions
		Ntrials(c) = M(c).ntrials;
	end
	% TotalTrials is total # of trials across conditions - this
	% will be used to allocate vectors for input into metricdist.mex function
	TotalTrials = sum(Ntrials);
	Times = cell(TotalTrials, 1);
	Labels = cell(TotalTrials, 1);
	Categories = zeros(TotalTrials, 1);

	tindex = 0;
	% loop through conditons
	for c = 1:Nconditions
		% loop through trials for this condition
		for n = 1:Ntrials(c)
			% assign spiketimes and labels (unit #) and category (condition)
			% to vectors used for metricdist.mex function
			tindex = tindex + 1;
			Times{tindex} = M(c).spikes{n};
			Labels{tindex} = uIndex*ones(size(M(c).spikes{n}));
			Categories(tindex) = c;
		end
	end

	% compute spike distance metric
	D = metricdist(1, Times, Labels, opts);

	% compute some stats from Ddistance (confusion) matrix
	[Dclust, H] = metric_cluster(D, CostVals, Categories, z);

	% store data
	Dlfh(uIndex).D = D;
	Dlfh(uIndex).Dclust = Dclust;
	Dlfh(uIndex).H = H;
end


% save to .mat file
save('spike_distance_output.mat', 'Dbbn', 'Dlfh', 'opts', 'z', '-MAT')


%{
% plot, for each unit, the information -vs- cost
for u = 1:length(Ddata)
	h = Ddata(u).Hvalue;
	figure(u)
	plot(CostVals, h, '.-')
	xlabel('cost/sec')
	ylabel('H (bits)')
	[a, aind] = max(h);
	hold on
		plot(CostVals(aind), a, 'ro')
	hold off
	title({	sprintf('File %s', filename1),	...
				sprintf('Unit %d', u + 1),			...
				sprintf('qmax=%.2f', CostVals(aind))	}, ...
				'Interpreter', 'none');
end

% store Ddata in OutData cell array
OutData{findex} = Ddata;
clear Ddata;
	


%}



