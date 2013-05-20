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

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
datapath = '/Users/sshanbhag/Work/Data/DataWave/batmat/AllSpkTrains_AllUnits';
stringfile = 'AllSpkTrains_AllUnits0to800_strings.mat';
syllfile = 'AllSpkTrains_AllUnits0to800_syllables.mat';

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = [0 2.^(-1:5)];
q = q(1:2);
% exponent
z = -2;
Nshuf = 50;		% # of shuffled/bootstrap tests

%-------------------------------------------------------------------------
% load spike trains 
%-------------------------------------------------------------------------
% MG settings:
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_testnewinfo_sec_0to100_strings.mat');
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_TestArray_2.mat');

% SJS 
load(fullfile(datapath, syllfile));

%-------------------------------------------------------------------------
% some values from data
%-------------------------------------------------------------------------
% # of trials, stimuli, units
% depends on test type
[nRows, Trials] = size(AllSpkTrains);
% # of conditions for data
% 13 for Type 1 (syllables) ; 8 for strings
Conditions = 13;
ConditionsToUse = 8;
nUnits = nRows/Conditions;
Ncost = length(q);

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% compare different methods for calculating spkd (spike distance metric)
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

% loop through units
for U = 1:1
	fprintf('********************************\n')
	fprintf('\tUNIT %d\n', U);
	fprintf('********************************\n')

	%-------------------------------------------------------------------------
	% extract the raw spiketrains for the current unit
	%-------------------------------------------------------------------------
	tmptrains = AllSpkTrains((U*Conditions+1 - Conditions):(U*Conditions),:);
	SpkTrains = tmptrains(1:ConditionsToUse, :);
	clear tmptrains;
	
	%-------------------------------------------------------------------------
	% settings for STAtoolkit
	% possible settings: {'plugin', 'tpmc', 'jack', 'ma', 'bub', 'chaoshen', 'ww',
	% 'nsb'}
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
	tic
	[STAresults{U}, STAshuf{U}] = metric_shuf(Sin, opts, Nshuf);
	% pull out value for different q values
	for Cindex = 1:Ncost
		fprintf('----------------------\n')
		fprintf('\t q = %.2f\n', opts.shift_cost(Cindex));
		fprintf('----------------------\n')	
		fprintf('\tH_sta: %.4f\n', STAresults{U}(Cindex).table.information.value);
		fprintf('\tConfusion Matrix:\n')
		% NOTE that confusion matrix from STA toolkit is TRANSPOSE of from from long
		% (slow) method!!!
		display_matrix(STAresults{U}(Cindex).cm', '%.1f');
	
		% get all the shuffled information bits
		hshuf = zeros(Nshuf, 1);
		for n = 1:Nshuf
			hshuf(n) = STAshuf{U}(Cindex, n).table.information.value;
		end		
		fprintf('\tShuffled Data:\n');
		fprintf('\t\tmean:\t%.4f\n', mean(hshuf));
		fprintf('\t\tstd:\t%.4f\n', std(hshuf));
		fprintf('\n')
	end
	fprintf('\n')
	sta_time = toc;
	fprintf('STAtoolkit method took %f seconds\n', sta_time);

end	% END U loop
