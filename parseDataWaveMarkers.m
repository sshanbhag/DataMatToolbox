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
	n
end

% store in Marker struct... might not be necessary, but keep it for now
Marker.M = M;

% Identify Stimulus types, store in Marker structure
Marker = identifyStimulus(Marker, dwinfo);

%-----------------------------------------------------------------------------
% Now check for different .wav file stimulus files
%-----------------------------------------------------------------------------
% Nunique			# of unique strings in Marker.text
% uniqueNames		unique values of strings in Marker.text, cell array
% uniqueIndices	indices where each of the unique strings are found in
%						Marker.WavFilename(L/R)
%
% if list of names is empty, leave wavFiles struct empty
%-----------------------------------------------------------------------------
% R channel
[names, indices, N] = findUniqueText(Marker.WavFilenameR);
% check if names is empty
if ~isempty(names)
	Marker.wavFilesR = struct(	'uniqueNames', names, ...
										'uniqueIndices', indices, ...
										'Nunique', N ...
										);
end
% L channel
[names, indices, N] = findUniqueText(Marker.WavFilenameL);
if ~isempty(names)
	Marker.wavFilesL = struct(	'uniqueNames', names, ...
										'uniqueIndices', indices, ...
										'Nunique', N ...
										);
end

%-----------------------------------------------------------------------------
% assign to output
%-----------------------------------------------------------------------------
out = Marker;
