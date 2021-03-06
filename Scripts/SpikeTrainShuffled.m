%------------------------------------------------------------------------
% SpikeTrainBatch
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed
%
%	uses code from expfun project by Chris Sumner as well
% 	as bits from the spike train toolbox
%------------------------------------------------------------------------
% Created: 1 December 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear old bits
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear all
% close all

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

% Cost for shifting spike, in units of 1/seconds
CostVals = [0 2.^(-2:9)];

% Time windows for spike trains.
SpikeWindowStart = 0;
SpikeWindowStop = 800;

%------------------------------------------------------------------------
% File Settings
%------------------------------------------------------------------------
dobjpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/dobj';
unitmatfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/UnitInfo.mat';
if ~exist('UnitInfo', 'var')
	load(unitmatfile);
end

return

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Read in File names
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if ~exist(fullfile(datapath, datalistfile), 'file')
	error('%s: file list %s not found', mfilename, fullfile(datapath, datalistfile));
else
	fprintf('\n\n%s: reading file list information from %s...\n\n\n', ...
					mfilename, fullfile(datapath, datalistfile));
end
fp = fopen(fullfile(datapath, datalistfile), 'r');
tmp = textscan(fp, '%s %s', 'Delimiter', ',');
fclose(fp);

% get # of files
Nfiles = length(tmp{1});

% build list of filenames
filelist = cell(Nfiles, 2);
for findex = 1:Nfiles
	filelist{findex, 1} = tmp{1}{findex};
	filelist{findex, 2} = tmp{2}{findex};
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

Ncost = length(CostVals);

% create list of error indices used to store problematic files
errorList = [];
nerrors = 0;

%------------------------------------------------------------------------
% loop through file list
%------------------------------------------------------------------------
for findex = 1:Nfiles
	% initialize ERRFLG that will be set to 1 if a file or data
	% error arises... this will allow processing for troublesome files to 
	% be skipped
	ERRFLG = 0;
	
	% get filenames for the long and short stim data files
	filename1 = filelist{findex, 1};
	filename2 = filelist{findex, 2};
	% create path+filename
	file1 = fullfile(datapath, filename1);
	file2 = fullfile(datapath, filename2);

	% read in data using load command (.efd files are mat files)
	data1 = load(file1, '-MAT');
	data2 = load(file2, '-MAT');

	% get spike train data - check in both raster and srcraster
	% fields of the data structures.  if data not found,
	% set ERRFLG to 1
	if isfield(data1, 'raster')
		data1Trains = data1.raster.ts_raster{1,1};
		TotalUnits1 = data1.raster.total_noof_units;	
	elseif isfield(data1, 'srcraster')
		data1Trains = data1.srcraster.ts_raster{1,1};
		TotalUnits1 = data1.srcraster.total_noof_units;	
	else
		warning('%s: could not find raster data for file %s', mfilename, file1);
		ERRFLG = 1;
	end
	if isfield(data2, 'raster')
		data2Trains = data2.raster.ts_raster{1,1};
		TotalUnits2 = data2.raster.total_noof_units;
	elseif isfield(data2, 'srcraster')
		data2Trains = data2.srcraster.ts_raster{1,1};
		TotalUnits2 = data2.srcraster.total_noof_units;
	else
		warning('%s: could not find raster data for file %s', mfilename, file2);
		ERRFLG = 1;
	end
	
	% check for matching numbers of units from the two files
	if TotalUnits1 ~= TotalUnits2
		% unit # mismatch - set ERRFLG
		warning('%s: # of units mismatch\n\t%s: %d units\n\t%s: %d units', mfilename, ...
					file1, TotalUnits1, file2, TotalUnits2);
		ERRFLG = 2;
	else
		TotalUnits = TotalUnits1;
	end
	
	% # of trials
	Trials1 = length(data1Trains(1,:));
	Trials2 = length(data2Trains(1,:));
	% make sure # of trials matches in the two files
	if Trials1 ~= Trials2
		warning('%s: # of trials mismatch\n\t%s: %d units\n\t%s: %d units', mfilename, ...
					file1, Trials1, file2, Trials2);
		ERRFLG = 3;
	else
		Trials = Trials1;
	end

	%------------------------------------------------------------------------
	% process data if no errors detected
	%------------------------------------------------------------------------
	if ~ERRFLG
		
		fprintf('Files: %s, %s\n', file1, file2);
		fprintf('# units: %d\n', TotalUnits);
		
		% Conditions (# of stimuli) = number of total conditions 
		% divided by the number of neurons	
		% get # of conditions
		Conditions1 = length(data1Trains(:,1))/TotalUnits1;
		Conditions2 = length(data2Trains(:,1))/TotalUnits2;

		% total # of conditions
		Conditions = Conditions1 + Conditions2;

		%Removes all spikes less than SpikeWindowStart and more than SpikeWindowStop
		for i = 1:Conditions1*TotalUnits1
			for i2 = 1:Trials1
				Index1 = find(data1Trains{i,i2} > SpikeWindowStart);
				data1Trains{i,i2} = data1Trains{i,i2}(Index1);       
				Index2 = find(data1Trains{i,i2} < SpikeWindowStop);
				data1Trains{i,i2} = data1Trains{i,i2}(Index2);
			end;    
		end;

		for i = 1:Conditions2*TotalUnits2
			for i2 = 1:Trials2
				Index1 = find(data2Trains{i,i2}>SpikeWindowStart);
				data2Trains{i,i2} = data2Trains{i,i2}(Index1);       
				Index2 = find(data2Trains{i,i2} < SpikeWindowStop);
				data2Trains{i,i2} = data2Trains{i,i2}(Index2);
			end;    
		end;

		% combine into AllSpkTrains
		nrows = Conditions1*TotalUnits1 + Conditions2*TotalUnits2;
		AllSpkTrains = cell(nrows, Trials);

		all_r = 0;
		for r = 1:Conditions1*TotalUnits1
			all_r = all_r + 1;
			for t = 1:Trials1
				AllSpkTrains{all_r, t} = data1Trains{r, t};
			end
		end
		for r = 1:Conditions2*TotalUnits2
			all_r = all_r + 1;
			for t = 1:Trials2
				AllSpkTrains{all_r, t} = data2Trains{r, t};
			end
		end

		
		
		
		
		% loop through units (start with 2nd unit, since unit(1) is garbage)
		unitIndex = 0;
		for Units = 2:TotalUnits
			fprintf('UNIT %d\n', Units);
			unitIndex = Units - 1;
			% loop through cost values
			for Cindex = 1:Ncost
				% get current cost value from list
				Cost = CostVals(Cindex);
				fprintf('\tCost = %f\n', Cost);

				% initialize Ddistance array
				Ddistance = zeros(Conditions,Conditions);
				
				% extracts the raw spiketrains for the units
				% to select a particular unit for example the first responsive one 
				% you can use units = 2
				SpkTrains = AllSpkTrains((Units*Conditions - Conditions+1):(Units*Conditions), :);
				

				% shuffle spike trains across Conditions
				for t = 1:Trials
					randindices(:, t) = randperm(Conditions);
				end
				for c = 1:Conditions
					randindices(c, :) = randindices(c, randperm(Trials));
				end
				SSpkTrains = SpkTrains(randindices);
				
				% loop for the number of conditions
				for C1 = 1:Conditions
					% loop for the number of trials;
					for T1 = 1:Trials

						% for each pair of conditions (stimuli), compute spike
						% distance metric for all trials
						for C2 = 1:Conditions
							for T2=1:Trials
								Dspike(C2,T2) = (spkd(SpkTrains{C1,T1}, SpkTrains{C2,T2}, Cost))^z;
							end;
						end;

						Dspike(C1,T1) = 0;   
						[a,b] = find(Dspike==Inf);
						Dspike(a,b) = 0;        
						Dsum = sum(Dspike,2);

						for C3 = 1:Conditions,
							if C3 == C1
								Dave(C3) = (Dsum(C3)/(Trials-1))^(1/z);
							else
								Dave(C3) = (Dsum(C3)/Trials)^(1/z);
							end;
						end;
						% find indices of minimum avg distance
						[c ,d] = min(Dave);
						Ddistance(C1,d) = Ddistance(C1,d) + 1;

						% clear some things
						clear Dspike Dsum Dave;
					end;
				end;

				% compute some stats from Ddistance (confusion) matrix
				[Hvalue, Hmean, Hprob] = Hstat(Ddistance);

				% ExsName= strcat(Filename.name,'-Cell-', num2str(Units))
				% xlswrite(ExsName, Ddistance)
				fprintf('\tHvalue: %f\n', Hvalue);
				fprintf('\tHmean: %f\n', Hmean);
				fprintf('\tHprob: %f\n', Hprob);
				fprintf('\tDdistance = \n');
				disp(Ddistance);
				
				% store data
				Ddata(unitIndex).Ddistance{Cindex} = Ddistance;
				Ddata(unitIndex).Hvalue(Cindex) = Hvalue;
				Ddata(unitIndex).Hmean(Cindex) = Hmean;
				Ddata(unitIndex).Hprob(Cindex) = Hprob;
				Ddata(unitIndex).UnitNumber = Units;
				Ddata(unitIndex).file1 = filelist{findex, 1};
				Ddata(unitIndex).file2 = filelist{findex, 2};
				
				clear Ddistance SpkTrains;
			end;
		end

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
	else
		% Error occurred
		fprintf('%s: ERROR %d DETECTED\n', mfilename, ERRFLG);
		fprintf('\t skipping to next files...\n\n')
		% increment error count and store error flag
		nerrors = nerrors + 1;
		errorList(1, 1) = findex;
		errorList(1, 2) = ERRFLG;
	end
	
end

% save to .mat file
save(outputmatfilename, '-MAT')




