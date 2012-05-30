function [out, errFlg] = identifyStimulus(Marker, dwinfo)
%------------------------------------------------------------------------
% [out, errFlg] = identifyStimulus(Markers, dwinfo)
%------------------------------------------------------------------------
% 
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
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 9 June, 2011 (SJS)
%
% Revisions:
%	30 May, 2012 (SJS): minor tweaks as Marker Object is built
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;
errFlg = 0;

Marker.StimulusTypeR = checkStimulusType(Marker, 'R');
Marker.StimulusTypeL = checkStimulusType(Marker, 'L');
out = Marker;
%------------------------------------------------------------------------


%------------------------------------------------------------------------
function StimulusType = checkStimulusType(Marker, C)
%------------------------------------------------------------------------
% StimulusType = checkStimulusType(Marker, C)
%------------------------------------------------------------------------
% check channel C stimulus type for data in Marker struct
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

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;

unknownFlag = 0;

% number of markers
N = length(Marker.(['WavFilename' C]));

% pre-allocate StimulusType cell array
StimulusType = cell(N, 1);

%----------------------------------------
% loop through number of markers
%----------------------------------------
for m = 1:N
	%----------------------------------------
	% first check if wav filename exists
	%----------------------------------------
	if isempty(Marker.(['WavFilename' C]){m})
		% no wav filename, set wavFlag to 0;
		wavFlag = 0;
	else
		wavFlag = 1;
	end
	
	%----------------------------------------
	% check for tone
	%----------------------------------------
if Marker.(['ToneFreq' C])(m) > 0
		toneFlag = 1;
	else
		toneFlag = 0;
	end
	
	%----------------------------------------
	% check for noise
	%----------------------------------------
	% low frequency
	if (Marker.(['BBNlowerFreq' C])(m) > 0)
		BBNlowFlag = 1;
	else
		BBNlowFlag = 0;
	end
	% high frequency
	if (Marker.(['BBNupperFreq' C])(m) > 0)
		BBNhiFlag = 1;
	else
		BBNhiFlag = 0;
	end
	
	if BBNlowFlag && BBNhiFlag
		% if both bbn low and bbn hi flags are set, broadband noise
		% stimulus was used; set noiseFlag to 1
% 		disp('noise stim detected')
		noiseFlag = 1;
	elseif (BBNlowFlag && ~BBNhiFlag) || (~BBNlowFlag && BBNhiFlag)
		% for some reason, only one of the noise freq. flags was set - 
		% deliver warning to user, set noiseFlag to 0;
		warning('%s: bizarre condition for stimulus m- either BBNlowerF or BBNupperF is zero', ...
						mfilename, m);
		unknownFlag = 1;
		noiseflag = 0;
	else
		noiseFlag = 0;
	end
	
	
	%----------------------------------------
	% now perform tests and determine stimulus type
	%----------------------------------------
	StimulusType{m} = [];

	% build vector from flags to simplify checking
	flagVector = [wavFlag toneFlag noiseFlag];

	if unknownFlag
		% stimulus type could not be determined
		StimulusType{m} = 'UNKNOWN';
		warning('%s: StimulusType%c{%d} is unknown', mfilename, C, m);

	elseif ~any(flagVector)
		% no sound for this channel
		StimulusType{m} = 'NO_SOUND';
	
	else
		% check if more than one flag is set
		if sum(flagVector) > 1 
			warning('%s: multiple stimuli for sweep # %d', m);
			disp(flagVector);
			multFlag = 1;
		else
			multFlag = 0;
		end
		
		if wavFlag
			StimulusType{m} = 'WAVFILE';
		end
		if toneFlag
			if multFlag
				StimulusType{m} = [StimulusType{m} 'TONE'];
			else
				StimulusType{m} = 'TONE';
			end
		end
		if noiseFlag
			if multFlag
				StimulusType{m} = [StimulusType{m} 'NOISE'];
			else
				StimulusType{m} = 'NOISE';
			end
		end
	end
			
end				
%------------------------------------------------------------------------
