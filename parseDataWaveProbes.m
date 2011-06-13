function [UnitData, ProbeData, errFlg] = parseDataWaveProbes(ProbeData, Marker)
%------------------------------------------------------------------------
% [UnitData, ProbeData, errFlg] = parseDataWaveProbes(ProbeData, Marker)
%------------------------------------------------------------------------
% parse datawave marker information
% 
% returns structure 
% 
%------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%
% 	errFlg	Error flag
% 					0		no error
% 					1		user cancelled file opening
% 					2		no fields found in header lines
% 					3		file not found
%
%------------------------------------------------------------------------
% See: readDataWaveHeader 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 3 June, 2011 (SJS)
% 	- uses code snipped from loadDWStimData.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: parse into epochs
%------------------------------------------------------------------------

errFlg = 0;

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

%-----------------------------------------------------------
% separate spike times into units
%-----------------------------------------------------------
% count number of probes
NProbes = length(ProbeData);

% initialize total unit counter
NUnits = 0;

% loop through probes
for p = 1:NProbes
	% find unique cluster id values
	unique_clusters = unique(ProbeData(p).cluster);
	
	% check if unique values were found
	if ~isempty(unique_clusters)
		% if so, figure out how many clusters there are
		ProbeData(p).nclusters = length(unique_clusters);
		
		% then, assign data to the UnitData() struct array that will hold
		% unit and timestamp data
		
		% initialize cluster index n
		n = 1;
		% loop through running count of total units, and store probe, unit number
		% indices and timestamps for this unit
		for u = (1 + NUnits) : (NUnits + ProbeData(p).nclusters)
			UnitData(u).probe = p;
			UnitData(u).unit = unique_clusters(n);
			UnitData(u).indices = find(ProbeData(p).cluster == UnitData(u).unit);
			UnitData(u).timestamp = ProbeData(p).t(UnitData(u).indices);
			n = n + 1;
		end
		
		NUnits = NUnits + ProbeData(p).nclusters;
	else
		ProbeData(p).nclusters = 0;
		error('%s: no units found for probe %d', mfilename, p);
	end
end


