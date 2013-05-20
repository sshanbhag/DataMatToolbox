%------------------------------------------------------------------------
% ViewSTAdata
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
% Created 16 April, 2013;
%
%------------------------------------------------------------------------

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
% build mat file name
[~, matfile] = fileparts(datafile);
matfile = [matfile '_SDMout.mat'];

%-------------------------------------------------------------------------
% load data
%-------------------------------------------------------------------------
load(matfile);

%-------------------------------------------------------------------------
% SDM settings
%-------------------------------------------------------------------------
Ncost = length(q);
Nunits = length(UnitList);

%% plot data
U = 1;
pFlag = 1;
while pFlag
	% plot H values
	plot(q, R(U).H, 'k.-')
	% find max, draw red circle at that point
	[maxval, maxind] = max(R(U).H);
	hold on
		plot(q(maxind), maxval, 'ro');
	hold off
	
	% plot shuffled mean, sd as area
	hold on
	plotsdarea(q, R(U).Hshuf_mean, R(U).Hshuf_sd)
	hold off
	
	legend({'H', 'Hmax', 'Hshuffled'})
	
	xlabel('q (Cost)')
	ylabel('Information (bits)')
	title({matfile, sprintf('Unit %d', UnitList(U))}, 'Interpreter', 'none')

	rFlag = 1;
	while rFlag
		resp = input('(z) previous, (x) next, unit #, (q) quit: ', 's');
		if ~isempty(resp)
			if ~isnan(str2double(resp))
				if (str2double(resp) > 0) && (str2double(resp) <= Nunits)
					U = str2double(resp);
					rFlag = 0;
				end
			elseif lower(resp) == 'z'
				if U > 1
					U = U - 1;
					rFlag = 0;
				end
			elseif lower(resp) == 'x'
				if U < Nunits
					U = U + 1;
					rFlag = 0;
				end
			elseif lower(resp) == 'q'
				rFlag = 0;
				pFlag = 0;
			end
		end
	end
end