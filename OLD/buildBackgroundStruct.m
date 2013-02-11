function Background = buildBackgroundStruct(Data)
%------------------------------------------------------------------------
% Background = buildBackgroundStruct(Data)
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
%
% Takes DataMat Data structure and splits up data by stimulus type
%
%------------------------------------------------------------------------
% Input Arguments:
%	Data				DataWave DATA structure (usu. from loadDWfile)
%
% Output Arguments:
%
% 	Background		structure of background (pre-first marker) data
% 	
%
%------------------------------------------------------------------------
% See also: buildStimulusStruct, loadDWfile, buildSpikes 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 4 August, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

% make local copy of Marker structure
Marker = Data.Marker;

%-----------------------------------------------------------------------------
% get the first marker time
%-----------------------------------------------------------------------------
initialMarkerTime = str2double(Marker.Timestamp(1));

if initialMarkerTime <= 0
	error('%s: initialMarkerTime (%f) <= 0 !!!!', mfilename, initialMarkerTime);
end

%------------------------------------------------------------------------
% Retrieve the background spikes for each unit
%------------------------------------------------------------------------
% make local copy of UnitData
UnitData = Data.UnitData;
Nunits = Data.Info.Nunits;

% allocate the Background structure
Background.Spiketimes = cell(Nunits, 1);
Background.initialMarkerTime = initialMarkerTime;
				
% loop through the units
for u = 1:Nunits
	% retrieve the spikes that are valid, using the initialMarkerTime as 
	% upper limit
	valid_times	= (UnitData(u).timestamp < initialMarkerTime);

	% store the values if spikes were found, 
	if any(valid_times)
		Background.Spiketimes{u} = UnitData(u).timestamp(valid_times);
	else
		%otherwise, set to empty array
		Background.Spiketimes{u} = [];
	end
end	% end of U loop



