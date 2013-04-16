%------------------------------------------------------------------------
% SpikeDistanceWithShuffle
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Used to compute Spike Distance Metric (Victor and Purpura) and
% a shuffled metric to establish confidence intervals
%------------------------------------------------------------------------
% 3 April, 2013, SJS
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear and cleanup
%------------------------------------------------------------------------
%------------------------------------------------------------------------
clear classes
close all

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% settings
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% power transformation for averaging, usually -2
z = -2;

% Cost of changing time
%CostVals = [0, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014,0.015, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5 ];
%CostVals=[0, 8, 9, 10, 11, 12, 13, 14, 15, 20, 30, 40, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500];
CostVals = [0, 2, 4, 8, 10, 20, 100];
Ncost = length(CostVals);

% # of shuffled data reps
Nbootstrap = 10;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Load spike train data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% MG settings:
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_testnewinfo_sec_0to100_strings.mat');
%load('E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2\AllSpkTrains_TestArray_2.mat');

% SJS 
load AllSpkTrains_test2

% # of trials, stimuli, units
% depends on test type
[nTrials, nReps]=size(AllSpkTrains);
nStimuli=8; %13 for Type 1 (syllables) ; 8 for strings
nUnits=nTrials/nStimuli;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% compute things
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% loop through units
%------------------------------------------------------------------------
for Units = 1%:nUnits
	% extracts the raw spiketrains for the current unit
	SpkTrains = AllSpkTrains((Units*nStimuli+1 - nStimuli):(Units*nStimuli),:);

	%------------------------------------------------------------------------
	% loop through cost values
	%------------------------------------------------------------------------
	for Cindex = 1:Ncost
        
		% get current cost value from list
		Cost = CostVals(Cindex);
		fprintf('Cost = %f\n', Cost);

		% initialize Ddistance array
		Ddistance= zeros(nStimuli,nStimuli);
      
		%------------------------------------------------------------------------
		% compute UNSHUFFLED metric
		%------------------------------------------------------------------------
		Ddistance = confusematrix(SpkTrains, Cost, z);		
		% compute some stats from Ddistance (confusion) matrix
		Hvalue = Hstat(Ddistance);

		% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
		% xlswrite(ExsName, Ddistance)
		fprintf('Unit %d\n', Units)
		fprintf('\tHvalue: %f\n', Hvalue);
		disp(Ddistance);
		Ddata(Units).Ddistance{Cindex} = Ddistance;
		Ddata(Units).Hvalue(Cindex) = Hvalue;

		%------------------------------------------------------------------------
		% compute SHUFFLED metric
		%------------------------------------------------------------------------
		fprintf('Computing bootstrap/shuffled values...');
		Hb = zeros(Nbootstrap, 1);
		for B = 1:Nbootstrap
			% shuffle spike trains across Stimuli and reps
			ShuffleTrains = shufflecell(SpkTrains);
			Dshuf = confusematrix(ShuffleTrains, Cost, z);
			% compute some stats from Ddistance (confusion) matrix
			Hb(B) = Hstat(Dshuf);
		end	% END B loop
		Ddata(Units).Hb{Cindex} = Hb;
		fprintf('... done!\n')
		fprintf('\t\tmean(Hb) = %f\n', mean(Hb));
		
	end	% END Cindex loop (cost values)
end	% END Nunits loop

return
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
    title({ sprintf('Unit %d', u),			...
            sprintf('qmax=%.2f', CostVals(aind))	}, ...
            'Interpreter', 'none');
end

% % store Ddata in OutData cell array
% OutData{findex} = Ddata;
% clear Ddata;
%objpath='E:\Bat - SingleCh Restrained\batmat\Converted Dobj files\ALL Dobj files2';
%save(fullfile(objpath,'TestFullDdata_sec_0to100_strings.mat'), 'Ddata', '-MAT');
