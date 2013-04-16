%------------------------------------------------------------------------
% SpikeDistanceSTA
%------------------------------------------------------------------------
% Spiketrain data are from 
%	stringfile = 'AllSpkTrains_AllUnits0to800_strings.mat';
%	syllfile = 'AllSpkTrains_AllUnits0to800_syllables.mat';
% 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% sharad shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created 15 April, 2013;
%
%------------------------------------------------------------------------



%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to toolbox
%------------------------------------------------------------------------
%------------------------------------------------------------------------
tmp = which('metric');
if isempty(tmp)
	fprintf('%s: setting up Spike paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/dev/Toolboxes/STAkitMAC'))
	fprintf('... done \n\n')
end
clear all;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
datapath = '/Users/sshanbhag/Work/Data/DataWave/batmat/AllSpkTrains_AllUnits';
stringfile = 'AllSpkTrains_AllUnits0to800_strings.mat';
syllfile = 'AllSpkTrains_AllUnits0to800_syllables.mat';

% select data input file
datafile = stringfile;
% build output file name
[~, outfile] = fileparts(datafile);
outfile = [outfile '_SDMout.mat'];

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
% cost (q) in 1/seconds
q = [0 2.^(-1:8)];
% exponent
z = -2;
UnitList = 0;
Nshuf = 100;		% # of shuffled tests
Ncost = length(q);

%-------------------------------------------------------------------------
% load spike trains 
%-------------------------------------------------------------------------
% MG settings:
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_testnewinfo_sec_0to100_strings.mat');
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_TestArray_2.mat');

% SJS 
load(fullfile(datapath, datafile));

%-------------------------------------------------------------------------
% some values from data
%-------------------------------------------------------------------------
% # of trials, stimuli, units
% depends on test type
[Nrows, Ntrials] = size(AllSpkTrains);
% # Conditions is 13 for syllables,  8 for strings
if strcmpi(datafile, syllfile)
	Conditions = 13;
else
	Conditions = 8; 
end
nTotalUnits = Nrows/Conditions;
% check for UnitList value out of range (or unspecified, <=0)
if isempty(UnitList) || (UnitList(1) <= 0)
	UnitList = 1:nTotalUnits;
	fprintf('Analyzing all %d units in file %s\n', nTotalUnits, datafile);
end
if any(UnitList > nTotalUnits)
	error('%s: UnitList value > nTotalUnits!', mfilename)
else
	Nunits = length(UnitList);
end

%-------------------------------------------------------------------------
% allocate storage for analysis
%-------------------------------------------------------------------------

% Store results in R struct
R = repmat( struct(	'H', zeros(Ncost, 1), ...			% H (information) values
							'CM', {cell(Ncost, 1)}, ...				% confusion matrices
							'Hshuf', zeros(Ncost, Nshuf), ...	% shuffled H vals
							'Hshuf_mean', zeros(Ncost, 1), ...	% avg Hshuf
							'Hshuf_sd', zeros(Ncost, 1) ...		% std. dev. Hshuf
						), ...
				Nunits, 1);
			
%-------------------------------------------------------------------------
% settings for STAtoolkit options
%-------------------------------------------------------------------------
% use the Treves-Panzeri- - correction for bias in entropy
% possible settings: {'plugin', 'tpmc', 'jack', 'ma', 'bub', 'chaoshen', 'ww',
% 'nsb'}
opts.entropy_estimation_method = {'tpmc'};
% shift_cost is set to values in q
opts.shift_cost = q;
% start and end time for analysis in each spike train - in seconds
opts.start_time = 0;
opts.end_time = 0.8;
% set metric_family to Dspike (0).  metric_family = 1 correspondes to Dinterval 
opts.metric_family = 0;
% set clustering exponent, z
opts.clustering_exponent = z;
% use parallel computation method, speeds up multiple q calculation
opts.parallel = 1;
% miscellaneous things
opts.unoccupied_bins_strategy = -1;
opts.possible_words = 'recommended';

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%% compare different methods for calculating spkd (spike distance metric)
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

fprintf('Using q values:\n');
fprintf('\t%.4f\n', q);
fprintf('\n');

U = 0;
for UnitIndex = UnitList
	% increment U index
	U = U + 1;
	fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
	fprintf('\tUNIT %d\n', UnitIndex);
	fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

	%-------------------------------------------------------------------------
	% extract the raw spiketrains for the current unit
	%-------------------------------------------------------------------------
	SpkTrains = AllSpkTrains(	(UnitIndex*Conditions+1 - Conditions) : ...
										(UnitIndex*Conditions), :);

	%-------------------------------------------------------------------------
	% then, use STA
	%-------------------------------------------------------------------------

	% build input using spikes2sta function
	Sin = spikes2sta(SpkTrains, 'timescale', 1, 'timeresolution', 1, ...
								'start_time', opts.start_time, 'end_time', opts.end_time);
	% call STAtoolkit function metric_shuf (same as metric, but also does shuffled
	% data procedure to check random floor
	[tmpr, tmps] = metric_shuf(Sin, opts, Nshuf);
	
	% pull out value for different q values
	for Cindex = 1:Ncost
		% store H value
		R(U).H(Cindex) = tmpr(Cindex).table.information.value;
		% store transposed confusion matrix -  
		% confusion matrix from STA toolkit is TRANSPOSE of that from long
		% (slow) method!!!
		R(U).CM{Cindex} = tmpr(Cindex).cm';
		% get all the shuffled information bits
		for n = 1:Nshuf
			R(U).Hshuf(Cindex, n) = tmps(Cindex, n).table.information.value;
		end
		% compute mean and std dev
		R(U).Hshuf_mean(Cindex) = mean(R(U).Hshuf(Cindex, :));
		R(U).Hshuf_sd(Cindex) = std(R(U).Hshuf(Cindex, :));
		% display information
		fprintf('----------------------\n')
		fprintf('\t q = %.2f\n', opts.shift_cost(Cindex));
		fprintf('----------------------\n')
		fprintf('\tConfusion Matrix:\n')
		for n = 1:Conditions
			fprintf('\t\t');
			fprintf('%.1f\t\t', R(U).CM{Cindex}(n, :));
			fprintf('\n');
		end
		fprintf('\tH (bits): %f\n', R(U).H(Cindex));
		fprintf('\tShuffled Data:\n');
		fprintf('\t\tmean:\t%f\n', R(U).Hshuf_mean(Cindex));
		fprintf('\t\tstd:\t%f\n', R(U).Hshuf_sd(Cindex));
	end
	fprintf('\n\n')

end	% END U loop

save(outfile, '-MAT', 'R', 'Conditions', 'UnitList', 'opts', ...
					'nTotalUnits', 'Nrows', 'Ntrials', 'Nshuf', 'q', 'z');


