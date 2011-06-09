function [out, errFlg] = parseDataWaveMarkers(M, dwinfo)
%------------------------------------------------------------------------
% [out, errFlg] = parseDataWaveMarkers(M)
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
% Created: 2 June, 2011 (SJS)
% 	- uses code snipped from loadDWStimData.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;
errFlg = 0;

%-----------------------------------------------------------------------------
% retrieve time stamps, values and text as arrays (or cell arrays)
%-----------------------------------------------------------------------------
%	Nmarkers								number of Datawave markers in .txt file
%	Marker.t(1:Nmarkers)				vector of Datawave marker timestamps
%	Marker.string{1:Nmarkers}		cell vector of marker strings
%	Marker.<tag name>(marker #)	numeric or cell vector that holds marker data
%-----------------------------------------------------------------------------

% check and store length of M struct array
Nmarkers = length(M);
if ~Nmarkers
	error('%s: no markers in marker struct array', mfilename);
end

% allocate arrays
Marker.string = cell(Nmarkers, 1);

for m = 1:MARKER_NMARKERS
	if any(strcmp(MARKER_TYPES{m}, {'int', 'float', 'double'}))
		Marker.(MARKER_TAGS{m}) = zeros(Nmarkers, 1);
	elseif strcmp(MARKER_TYPES{m}, 'char')
		Marker.(MARKER_TAGS{m}) = cell(Nmarkers, 1);
	else
		error('%s: undefined marker type %s', mfilename, MARKER_TYPES{n});
	end
end

% loop through markers (in M() struct array), pulling out text and value
for n = 1:Nmarkers
	for m = 1:MARKER_NMARKERS
		if any(strcmp(MARKER_TYPES{m}, {'int', 'float', 'double'}))
			Marker.(MARKER_TAGS{m})(n) = str2num(M(n).string{m});
		elseif strcmp(MARKER_TYPES{m}, 'char')
			Marker.(MARKER_TAGS{m}){n} = M(n).string{m};
		else
			error('%s: undefined marker type %s', mfilename, MARKER_TYPES{n});
		end		
	end
	Marker.string{n} = M(n).string;
end

Marker.M = M;

% save markertest.mat Marker dwinfo -MAT
% 
% keyboard

%-----------------------------------------------------------------------------
% Now check for different .wav file stimulus files
%-----------------------------------------------------------------------------
%	algorithm: 
% 		1) obtain stimulus file name (searchstr) from list of 
% 			testvals
% 		2) find locations of searchstr in full list Marker.text
% 		3) see if this string already exists in the uniqueText list
% 				If yes:
% 					increment the testIndex so that the next value in
% 					testvals becomes the new searchstr
% 				If no:
% 					increment NuniqueStim count
% 					store the current searchstr in uniqueText
% 					store indices of searchstr in Marker.text in the uniqueIndices list
% 					rebuild the testvals list by removing instances of the current
% 						searchstr
% 						
% *Will need a different (slightly) approach for frequency/response
% area data
%-----------------------------------------------------------------------------
% NuniqueStim		# of unique strings in Marker.text
% uniqueText		unique values of strings in Marker.text, cell array
% uniqueIndices	indices where each of the unique strings are found in Marker.text
% testvals			temporary storage of unmatched strings
%-----------------------------------------------------------------------------

% Initialize some arrays and variables
NuniqueStim = 0;
uniqueText = {};
uniqueIndices = {};
textIndex = 1;
testIndex = 1;
runFlag = 1;


% Make a copy of the Marker.text cell array, stored in testvals to search
testvals = Marker.text;

% loop through markers or until runFlag is set to 0 by check within the loop
while runFlag && (NuniqueStim <= Nmarkers)
	% grab new string value to search for
	searchstr = testvals{testIndex};
	% look for matches in original Marker.text array
	matchindices = strcmp(searchstr, Marker.text);
	
	% check if any matches were found
	if ~sum(matchindices)
		% no matches found 
		disp(['no match found for ' searchstr])
		% increment testIndex
		testIndex = testIndex + 1;
		% check if testIndex is out of bounds
		if testIndex > length(testvals)
			% stop loop
			runFlag = 0;
		end
		
	else
		% match of searchstr found in testvals{}

		% check if new value found by searching for string in uniqueText
		% cell array (which is where the unique strings are stored)
		uniqueCheck = find(strcmp(searchstr, uniqueText));
		
		if ~isempty(uniqueCheck)
			% string is not new, so skip ahead to a new test value
			testIndex = testIndex + 1;
			% if new testIndex value is out of bounds, stop loop
			if testIndex > length(testvals)
				runFlag = 0;
			end
		else
			% New (unique) value found
			
			% increment NuniqueStim counter
			NuniqueStim = NuniqueStim + 1;
			% store the current search string in uniqueText cell array
			uniqueText{NuniqueStim, 1} = searchstr;
			% save Marker.text indices for this search string
			uniqueIndices{NuniqueStim, 1} = find(matchindices);
		
			% rebuild the testvals cell array by removing the instances of the
			% current search string
			newIndices = find(~matchindices);
			if isempty(newIndices)
				% if no new strings were found, stop search
				runFlag = 0
			else
				% otherwise, assign unique strings to testvals
				testvals = cell(length(newIndices), 1);
				for n = 1:length(newIndices)
					testvals{n} = Marker.text{newIndices(n)};
				end
				testIndex = 1;
			end
			
		end	% end of "if ~isempty(uniqueCheck)"
	end	% end of "if ~sum(matchindices)"
end	% end of "while runFlag && (NuniqueStim <= Nmarkers)"


















out = Marker;

