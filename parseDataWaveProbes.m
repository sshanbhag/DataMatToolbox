function [UnitData, ProbeData, errFlg] = parseDataWaveProbes(ProbeData, Marker)
%------------------------------------------------------------------------
% [UnitData, ProbeData, errFlg] = parseDataWaveProbes(ProbeData, Marker)
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
% parse datawave marker information
% 
%------------------------------------------------------------------------
% Input Arguments:
%	ProbeData			structure containing spike information
% 	Marker				structure containing event information
% 
% Output Arguments:
% 	UnitData(n)
% 		probe				probe ID #
% 		unit				unit ID # (usually equal to n-1)
% 		indices			vector of indexes into master event list [1xN double]
% 		timestamp		vector of time stamps, microseconds [1xN double]
% 		sorted			1 if spikes are sorted/thresholded, 0 for unsorted spikes
%
% 	ProbeData	spike information structure, cluster information added
% 		t					timestamp vector (microseconds)
% 		cluster			cluster ID number for spikes 
% 								0 indicates unsorted spike
% 		nclusters		number of clusters on this probe (electrode or tetrode)
% 	
%	errFlg		Error flag
% 		0					no error
% 		1					user cancelled file opening
% 		2					no fields found in header lines
% 		3					file not found
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
% 	5 July, 2011 (SJS)
% 	 -	added sorted parameter to UnitData struct
% 	 -	added comments
%------------------------------------------------------------------------
% TO DO: parse into epochs
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

errFlg = 0;

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
			if UnitData(u).unit > 0
				UnitData(u).sorted = 1;
			else
				UnitData(u).sorted = 0;
			end
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


