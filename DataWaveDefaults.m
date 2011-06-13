%------------------------------------------------------------------------
% DataWaveDefaults
%------------------------------------------------------------------------
% script that defines default values for DataWave Text files
% 
%------------------------------------------------------------------------
% See: readDataWaveHeader 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 25 May, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


%-----------------------------------------------------------
% default # of header lines
%-----------------------------------------------------------
N_HEADER_LINES = 1;


% Marker format information

MARKER_TAGS = { ...
	'Timestamp', ...
	'id', ...
	'OutputFilename', ...
	'AttenuationR', ...
	'BBNlowerFreqR', ...
	'BBNupperFreqR', ...
	'AmplitudeR', ...
	'TimeShiftR', ...
	'RampUpR', ...
	'HoldTimeR', ...
	'RampDownR', ...
	'OutputTimestampR', ...
	'OutputTimeWithDelayR', ...
	'FixedDelayR', ...
	'PA5idR', ...
	'WavFilenameR', ...
	'ToneFreqR', ...
	'PhaseDegR', ...
	'AttenuationL', ...
	'BBNlowerFreqL', ...
	'BBNupperFreqL', ...
	'AmplitudeL', ...
	'TimeShiftL', ...
	'RampUpL', ...
	'HoldTimeL', ...
	'RampDownL', ...
	'OutputTimestampL', ...
	'OutputTimeWithDelayL', ...
	'FixedDelayL', ...
	'PA5idL', ...
	'WavFilenameL', ...
	'ToneFreqL', ...
	'PhaseDegL' ...
};

MARKER_TYPES = { ...
	'char', ...
	'int', ...
	'char', ...
	'float', ...
	'float', ...
	'float', ...
	'float', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'char', ...
	'float', ...
	'float', ...
	'float', ...
	'float', ...
	'float', ...
	'float', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'int', ...
	'char', ...
	'float', ...
	'float', ...
};

MARKER_UNITS = { ...
	'usec', ...
	'', ...
	'', ...
	'dB', ...
	'Hz', ...
	'Hz', ...
	'mV', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'', ...
	'', ...
	'Hz', ...
	'Degrees', ...
	'dB', ...
	'Hz', ...
	'Hz', ...
	'mV', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'usec', ...
	'', ...
	'', ...
	'Hz', ...
	'Degrees' ...
};

MARKER_INFO = [MARKER_TAGS' MARKER_TYPES' MARKER_UNITS'];

MARKER_NMARKERS = length(MARKER_TAGS);

STIMULUS_TYPES = { ...
	'NOISE', ...
	'WAVFILE', ...
	'TONE' ...
	'NO_SOUND', ...
	'UNKNOWN' ...
};

TEST_TYPES = { ...
	'RATE_LEVEL', ...
	'FREQUENCY_RESPONSE', ...
	'NOISE_REPONSE' ...
};




