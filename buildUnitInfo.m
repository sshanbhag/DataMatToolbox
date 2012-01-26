function UnitInfo = buildUnitInfo(D, Stimulus)
%------------------------------------------------------------------------
% buildUnitInfo
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
% builds a struct array that contains information about the probe and
% cluster # from the spikesorting algorithm
% 
%------------------------------------------------------------------------
% Input Arguments:
%	D		Datawave struct
%	S		Stimulus struct
% 
% Output Arguments:
%	UnitInfo		[# units, 1] unit information struct array with fields:
% 		probe			Probe # for this unit
% 		cluster		cluster # for this unit (ranges from 0-n)
% 		unit			unit #
%
%------------------------------------------------------------------------
% See:  
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 	25 January 2012 (SJS):
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%-----------------------------------------------------------
% initialize some things
%-----------------------------------------------------------
nProbes = length(D.Probe);
tmp = size(Stimulus(1).Spiketimes);
nUnitsFromStimulus = tmp(1);
nUnits = 0;

%-----------------------------------------------------------
% now build the UnitInfo structure array
%-----------------------------------------------------------

for p = 1:nProbes
	if D.Probe(p).nclusters > 0
		for c = 1:D.Probe(p).nclusters
			nUnits = nUnits + 1;
			UnitInfo(nUnits).probe = p;
			UnitInfo(nUnits).cluster = c - 1;
			UnitInfo(nUnits).unit = nUnits;
		end
	end
end

if nUnits == 0
	warning('%s: no units found!', mfilename);
	UnitInfo = [];
end



