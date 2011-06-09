% function [out, errFlg] = identifyStimulus(Markers, dwinfo)
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
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 9 June, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


load markertest.mat

%-----------------------------------------------------------
% load defaults
%-----------------------------------------------------------
DataWaveDefaults;
errFlg = 0;


% check R channel stimulus type
for m = 1:length(Marker.WavFilenameR)
	% first check if wav filename exists
	if isempty(Marker.WavFilenameR{m})
		% no wav filename, so check other possibilities
		
		% check for tone
		if Marker.ToneFreqR(m) > 0
			Marker.StimulusTypeR{m} = 'TONE';
		elseif Marker.BBNlowerFreqR(m) > 0
			Marker.StimulusTypeR{m} = 'NOISE';
		else
			Marker.StimulusTypeR{m} = 'UNKNOWN';
		end
	else
		Marker.StimulusTypeR{m} = 'WAVFILE';
	end
end

% check L channel stimulus type
for m = 1:length(Marker.WavFilenameL)
	% first check if wav filename exists
	if isempty(Marker.WavFilenameL{m})
		% no wav filename, so check other possibilities
		
		% check for tone
		if Marker.ToneFreqL(m) > 0
			Marker.StimulusTypeL{m} = 'TONE';
		elseif Marker.BBNlowerFreqL(m) > 0
			Marker.StimulusTypeL{m} = 'NOISE';
		else
			Marker.StimulusTypeL{m} = 'UNKNOWN';
		end
	else
		Marker.StimulusTypeL{m} = 'WAVFILE';
	end
end

