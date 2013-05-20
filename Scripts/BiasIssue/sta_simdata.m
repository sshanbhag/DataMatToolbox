%------------------------------------------------------------------------
% scratch
%------------------------------------------------------------------------
% examining issues with STAtoolkit and shuffled data
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% sharad shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created 13 April, 2013;
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to toolboxes
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% STAkitMAC
tmp = which('metric');
if isempty(tmp)
	fprintf('%s: setting up STAtoolkit paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/STAkitMAC'))
	fprintf('... done \n\n')
end
clear tmp;
% SDM tools (with spkdl MEX file)
tmp = which('spkdl');
if isempty(tmp)
	fprintf('%s: setting up SDM paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/SDM'))
	fprintf('... done \n\n')
end
clear tmp;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% misc things
%------------------------------------------------------------------------
%------------------------------------------------------------------------
sepstr = '---------------------------------------------------------------';

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = [0 2.^(-1:10)];
q = q(1:4);
% exponent
z = -2;
Nshuf = 50;		% # of shuffled/bootstrap tests

%-------------------------------------------------------------------------
% Simulation values
%-------------------------------------------------------------------------
% # trials per condition
Trials = 30;
% # of conditions for data
Conditions = 13;
% # of cost vals
Ncost = length(q);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% create simulated data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% SpkTrains = build_simdata(	'Categories', Conditions, 'Reps', Trials, ...
% 									'Nspikes', randi([0 10], Conditions, 1), ...
% 									'Random_ISI', 10);

%------------------------------------------------------------------------
% Spike trains with following characteristics:
%	Mean # spikes varies with Condition
%	MeanISI is constant for all conditions.
%------------------------------------------------------------------------
% allocate SpkTrains cell matrix
SpkTrains1 = cell(Conditions, Trials);
% set mean # spikes per condition - this will be used to generate a random
% number of spikes using either gaussian or poisson distribution
SpikeMeanPerCondition =  (1:Conditions);
% Mean ISI per condition (in milliseconds)
MeanISIPerCondition = 10 * ones(Conditions, 1);
% pre-allocate Nspikes for all conditions and trials
Nspikes = zeros(Conditions, Trials);
% loop through conditions
for c = 1:Conditions
	% generate random number of spikes for each trial (using spikemean)
	Nspikes(c, :) = poissrnd(SpikeMeanPerCondition(c), Trials, 1);
	% loop through trials
	for t = 1:Trials
		SpkTrains{c, t} = poisson_spiketrain(Nspikes(c, t), MeanISIPerCondition(c));
	end	% END Trials loop
end	% END Conditions loop


%------------------------------------------------------------------------
% Spike trains with following characteristics:
%	# spikes is random across conditions, in interval [0 10]
%	MeanISI is constant for all conditions.
%------------------------------------------------------------------------
% allocate SpkTrains cell matrix
SpkTrains2 = cell(Conditions, Trials);
% Mean ISI per condition (in milliseconds)
MeanISIPerCondition2 = 10 * ones(Conditions, 1);
% generate random number of spikes for each trial
Nspikes2 = randi([0 10], Conditions, Trials);
% loop through conditions
for c = 1:Conditions
	% loop through trials
	for t = 1:Trials
		SpkTrains2{c, t} = poisson_spiketrain(Nspikes2(c, t), MeanISIPerCondition2(c));
	end	% END Trials loop
end	% END Conditions loop


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
% run analysis
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


SpkTrains = SpkTrains2;

%-------------------------------------------------------------------------
% settings for STAtoolkit
% possible settings: {'plugin', 'tpmc', 'jack', 'ma', 'bub', 'chaoshen', 'ww',
% 'nsb'}
%-------------------------------------------------------------------------
opts.entropy_estimation_method = {'plugin'};
opts.shift_cost = q;
opts.start_time = 0;
opts.end_time = .8;
opts.metric_family = 0;
opts.clustering_exponent = z;	
opts.parallel = 1;
opts.unoccupied_bins_strategy = -1;
opts.possible_words = 'recommended';
Ncost = length(opts.shift_cost);

% build input using spikes2sta function
Sin = spikes2sta(SpkTrains, 'timescale', 1, 'timeresolution', 1, ...
							'start_time', opts.start_time, 'end_time', opts.end_time);						
						
% call STAtoolkit function metric_shuf (same as metric, but also does shuffled
% data procedure to check random floor
[STAresults, STAshuf] = metric_shuf(Sin, opts, Nshuf);

% pull out value for different q values

% allocate storage
data.q = q;
data.CM = cell(Ncost, 1);
data.Hsta = zeros(Ncost, 1);
data.Hshuf.values = cell(Ncost, 1);
data.Hshuf.mean = zeros(Ncost, 1);
data.Hshuf.sd = zeros(Ncost, 1);

for Cindex = 1:Ncost
	fprintf('----------------------\n')
	fprintf('\t q = %.2f\n', opts.shift_cost(Cindex));
	fprintf('----------------------\n')	
	data.CM{Cindex} = STAresults(Cindex).cm';
	data.Hsta(Cindex) = STAresults(Cindex).table.information.value;
	fprintf('\tH_sta: %.4f\n', data.Hsta(Cindex));
	fprintf('\tConfusion Matrix:\n')
	% NOTE that confusion matrix from STA toolkit is TRANSPOSE of from from long
	% (slow) method!!!
	display_matrix(data.CM{Cindex}, '%.1f');

	% get all the shuffled information bits
	data.Hshuf.values{Cindex} = zeros(Nshuf, 1);
	for n = 1:Nshuf
		data.Hshuf.values{Cindex}(n) = STAshuf(Cindex, n).table.information.value;
	end
	data.Hshuf.mean(Cindex) = mean(data.Hshuf.values{Cindex});
	data.Hshuf.sd(Cindex) = std(data.Hshuf.values{Cindex});
	fprintf('\tShuffled Data:\n');
	fprintf('\t\tmean:\t%.4f\n', data.Hshuf.mean(Cindex));
	fprintf('\t\tstd:\t%.4f\n', data.Hshuf.sd(Cindex));
	fprintf('\n')
end
fprintf('\n')


plot(q, data.Hsta, 'r.-');
hold on
plotsdarea(q, data.Hshuf.mean, data.Hshuf.sd);
hold off