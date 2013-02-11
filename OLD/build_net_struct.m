function Sout = build_net_struct(Sdata, Slist)
%-----------------------------------------------------------------------------
% Sout = build_net_struct(Sdata, Slist)
%-----------------------------------------------------------------------------
% 
% Description
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	Input		input info
% 
% Output Arguments:
% 	Output	output info
%
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 8 February, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% check inputs
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% get # of files, 
%-----------------------------------------------------------------------------
Nlist = length(Slist);
% get the # of conditions from the length of the data indices in Slist{1, 1}
Nconditions = length(Slist{1, 1});

%-----------------------------------------------------------------------------
% loop through list of valid data (Slist)
%-----------------------------------------------------------------------------
% initialize sindx counter to 0
sindx = 0;

for f = 1:Nlist
	
	% get needed indices from the valid data list
	validDataIndices = Slist{f, 1};
	validUnitIndices = Slist{f, 5};
	validAttenIndices = Slist{f, 7};
	[NvalidUnits, tmp] = size(validUnitIndices);
	
	%-------------------------------------------------------------------------
	% find indices to non-zero clusters (i.e., sorted clusters only)
	%-------------------------------------------------------------------------
	clustn = zeros(NvalidUnits, Nconditions);
	for u = 1:NvalidUnits
		for c = 1:Nconditions
			clustn(u, c) = Sdata{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
		end
	end
	nonZeroUnits = find(clustn(:, 1));
	Nunits = length(nonZeroUnits);
		
	%-------------------------------------------------------------------------
	% loop through units
	%-------------------------------------------------------------------------
	for u = 1:Nunits
		% get the index value for this non-zero unit from the nonZeroUnits list
		uIndex = nonZeroUnits(u);

		% loop through conditions
		lfFlag = 0;
		for c = 1:Nconditions
			
			% get Datum for this condition
			D = Sdata{validDataIndices(c)};
			% columns in validUnitIndices correspond to conditions/files
			U = D.UnitData(validUnitIndices(u, c));
			attIndex = validAttenIndices(1, c);
			
			% check to make sure that the cluster # ~= 0
			if U.UnitInfo.cluster ~= 0
				% if so, store the data
				
				% increment output struct index
				sindx = sindx + 1;
				%-----------------------------------------------------------------
				% file				original data file name
				% condition			stimulus condition (odor)
				% unit				unit # (from loadDWfile...)
				% probe				probe #
				% cluster			cluster #
				% BG_timestamps	timestamps for background activity for this unit
				%						(from first 5 seconds of continuous data collected)
				% atten				stimulus attenuation (dB)
				% ntrials			# of trials for this stimulus (a.k.a. # sweeps)
				% spikes				cell array of neural spike times in milliseconds
				%-----------------------------------------------------------------
				Sout(sindx).file			= D.Info.file;
				Sout(sindx).condition	= D.Info.condition;
				Sout(sindx).unit			= U.UnitInfo.unit;
				Sout(sindx).probe			= U.UnitInfo.probe;
				Sout(sindx).cluster		= U.UnitInfo.cluster;
				Sout(sindx).BG_tstamps	= U.BG_timestamps;
				Sout(sindx).atten			= D.AttenVals(attIndex);
				Sout(sindx).ntrials		= U.ntrials(attIndex);
				Sout(sindx).spikes		= U.spikes{attIndex}; 

			end	% end if
		end		% END c LOOP
	end		% END u LOOP
end		%END f LOOP


