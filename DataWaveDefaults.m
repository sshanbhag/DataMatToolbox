function DataWaveDefaults(varargin)
%------------------------------------------------------------------------
% DataWaveDefaults
%------------------------------------------------------------------------
% DataMat toolbox
%------------------------------------------------------------------------
% script that defines default values for DataWave Text files
% 
% Values defined:
% 
% 	MARKER_TAGS			are the Event Data values from the DDF files
% 							as seen by the NeuroShare API
% 		Values:
% 			SoundType<R/L>
% 			Attenuation<R/L>
% 			WavFilename<R/L>
% 			BBNlowerFreq<R/L>
% 			BBNupperFreq<R/L>
% 			Amplitude<R/L>
% 			TimeShift<R/L>
% 			RampUp<R/L>
% 			HoldTime<R/L>
% 			RampDown<R/L>
% 			OutputTimestamp<R/L>
% 			OutputTimeWithDelay<R/L>
% 			FixedDelay<R/L>
% 			PA5id<R/L>
% 			ToneFreq<R/L>
% 			PhaseDeg<R/L>
% 			OutputFile
% 
% 	MARKER_TYPES		specifies the type of marker information -
% 		Values:
% 			'char'			character or string
% 			'int'				integer
% 			'float'			floating point
% 
% 	MARKER_UNITS		specifies units of marker information.
% 		Values:
% 			'usec'			microseconds
% 			'dB'				decibels
% 			'Hz'				Hertz
% 			'mV'				millivolts
% 			'Degrees'		degrees (angle)
% 			''					blank/dimensionless
% 
% 	MARKER_NMARKERS	is # of marker elements 
% 							(number of MARKER_TAGS, TYPES and UNITS)
% 
% 	MARKER_INFO			is a [MARKER_NMARKERS X 3] cell array of values with columns
% 							corresponding to MARKER_TAGS, MARKER_TYPES, and MARKER_UNITS
% 							i.e., MARKER_INFO = [MARKER_TAGS' MARKER_TYPES' MARKER_UNITS']
%------------------------------------------------------------------------
% See: DWdata, NS objects 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 25 May, 2011 (SJS)
%
% Revisions:
% 	8 August, 2011 (SJS):
%	 -	added pre/post time
%	30 May, 2012 (SJS):
% 	 -	updated comments/docs
% 	 - functionalized
%	14 Jan 2013 (SJS)
%	 - updated for use with NeuroShare event markers
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%-----------------------------------------------------------
% figure out what to do
%-----------------------------------------------------------
% see if mat file exists; if so, load that (in caller's
% workspace) and return
if isempty(varargin)
	if exist('DataWaveDefaults.mat', 'file')
		evalin('caller', 'load(''DataWaveDefaults.mat'')')
		return
	end
end
% otherwise, go ahead and create the default vars.
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
% default plot options
%-----------------------------------------------------------
% time before 
PRE_PLOTTIME = 250;
% and after sweep to plot 
POST_PLOTTIME = 250;
%-----------------------------------------------------------
%-----------------------------------------------------------
N_CHANNELS = 2;
L = 1;
R = 2;

%-----------------------------------------------------------
%-----------------------------------------------------------
% Default Marker Base values (minus L/R designation
%-----------------------------------------------------------
% Column 1: Marker Base tags
% Column 2: Marker data type
% Column 3: Marker data units
% 
% Marker Base Tags are the Marker Event CSV tag without the 
% 'R' and 'L' channel character.
% Marker data types:
% 	'char'		character or string
% 	'int'			integer
% 	'float'		floating point
% Marker data units specifies units of marker information.
%	'usec'			microseconds
% 	'dB'				decibels
% 	'Hz'				Hertz
% 	'mV'				millivolts
% 	'Degrees'		degrees (angle)
% 	''					blank/dimensionless	
%-----------------------------------------------------------
%-----------------------------------------------------------
MARKER_BASE =	{...
	'SoundType'					'char'	''				;	...
	'Attenuation'				'float'	'dB'			;	...
	'WavFilename'				'char'	''				;	...
	'BBNlowerFreq'				'int'		'Hz'			;	...
	'BBNupperFreq'				'int'		'Hz'			;	...
	'Amplitude'					'int'		'mV'			;	...
	'TimeShift'					'int'		'msec'		;	...
	'RampUp'						'int'		'msec'		;	...
	'HoldTime'					'int'		'msec'		;	...
	'RampDown'					'int'		'msec'		;	...
	'OutputTimestamp'			'int'		'usec'		;	...
	'OutputTimeWithDelay'	'int'		'usec'		;	...
	'FixedDelay'				'int'		'msec'		;	...
	'PA5id'						'int'		''				;	...
	'ToneFreq'					'int'		'Hz'			;	...
	'PhaseDeg'					'int'		'degrees'	;	...
	'OutputFile'				'char'	''				;	...
};
%-----------------------------------------------------------
%-----------------------------------------------------------
% Marker format information
%-----------------------------------------------------------
%-----------------------------------------------------------
% MARKER_TAGS are the field names for the Event Markers
[MARKER_NBASE, ncols] = size(MARKER_BASE);
MARKER_TAGS = cell(2*MARKER_NBASE, 1);
for n = 1:MARKER_NBASE
	% assign tags for right channel
	MARKER_TAGS{n} = [MARKER_BASE{n, 1} 'R'];
	% assign tags for left channel
	% !!! note that left channel is presumed to NOT have
	% !!! an 'OutputFile' tag/field... downstream funcitons
	% !!! should take this into account!!!!!!!
	MARKER_TAGS{n+MARKER_NBASE} = [MARKER_BASE{n, 1} 'L'];
end
% MARKER_TYPES specifies the type of marker information
MARKER_TYPES = [MARKER_BASE(:, 2); MARKER_BASE(:, 2)];
% MARKER_UNITS specifies units of marker information.
MARKER_UNITS = [MARKER_BASE(:, 3); MARKER_BASE(:, 3)];
% MARKER_NMARKERS = length (number) of MARKER_TAGS, TYPES and UNITS
MARKER_NMARKERS = length(MARKER_TAGS);
% MARKER_INFO is a [MARKER_NMARKERS X 3] cell array of values with columns 
% corresponding to MARKER_TAGS, MARKER_TYPES, and MARKER_UNITS
MARKER_INFO = [MARKER_TAGS MARKER_TYPES MARKER_UNITS];
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
STIMULUS_TYPES = { ...
	'NOISE', ...
	'WAVFILE', ...
	'TONE' ...
	'NO_SOUND', ...
	'BACKGROUND', ...
	'UNKNOWN' ...	
};
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
TEST_TYPES = { ...
	'RATE_LEVEL', ...
	'FREQUENCY_RESPONSE', ...
	'NOISE_REPONSE', ...
	'BACKGROUND' ...
};
%-----------------------------------------------------------
%-----------------------------------------------------------


%-----------------------------------------------------------
%-----------------------------------------------------------
STIMULUS_TAGS = { ...
	'StimulusType', ...
	'Attenuation', ...
	'BBNlowerFreq', ...
	'BBNupperFreq', ...
	'Amplitude', ...
	'TimeShift', ...
	'RampUp', ...
	'HoldTime', ...
	'RampDown', ...
	'WavFilename', ...
	'ToneFreq', ...
	'PhaseDeg', ...
};
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
STIMULUS_SEARCH_TAGS = { ...
	'StimulusType', ...
	'BBNlowerFreq', ...
	'BBNupperFreq', ...
	'Amplitude', ...
	'TimeShift', ...
	'RampUp', ...
	'HoldTime', ...
	'RampDown', ...
	'WavFilename', ...
	'ToneFreq', ...
	'PhaseDeg', ...
};
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
% get number of tags used to identify stimuli that are defined in the
% STIMULUS_SEARCH_TAGS{} vector
NTAGS_PER_CHANNEL = length(STIMULUS_SEARCH_TAGS);
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
TONE_VAR_TAGS = { ...
 	'Attenuation', ...
	'ToneFreq', ...
	'PhaseDeg' ...
};
NOISE_VAR_TAGS = { ...
 	'Attenuation', ...
	'BBNlowerFreq', ...
	'BBNupperFreq' ...	
};
WAV_VAR_TAGS = { ...
 	'Attenuation', ...
	'WavFilename' ...	
};
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
%-----------------------------------------------------------
STIMULUS_STRUCT_FIELDS = {	...
	'Type',			...
	'Channel',		...
	'Indices',		...
	'Var'
};

% merge in the marker tags
STIMULUS_STRUCT_FIELDS = [STIMULUS_STRUCT_FIELDS MARKER_TAGS'];
%-----------------------------------------------------------
%-----------------------------------------------------------

%-----------------------------------------------------------
% save the variables in mat file
%-----------------------------------------------------------
fprintf('Saving file DataWaveDefaults.mat...')
save('DataWaveDefaults.mat', ...
	'N_CHANNELS', ...
	'PRE_PLOTTIME', ...
	'POST_PLOTTIME', ...
	'MARKER_BASE', ...
	'MARKER_NBASE', ...
	'MARKER_TAGS', ...
	'MARKER_TYPES', ...
	'MARKER_UNITS', ...
	'MARKER_INFO', ...
	'MARKER_NMARKERS', ...
	'STIMULUS_TYPES', ...
	'TEST_TYPES', ...
	'STIMULUS_TAGS', ...
	'STIMULUS_SEARCH_TAGS', ...
	'NTAGS_PER_CHANNEL', ...
	'TONE_VAR_TAGS', ...
	'NOISE_VAR_TAGS', ...
	'WAV_VAR_TAGS', ...
	'STIMULUS_STRUCT_FIELDS', ...
	'L', ...
	'R', ...
	'-MAT' ...
);

% and load into caller's workspace
fprintf('...and loading variables\n');

evalin('caller', 'load(''DataWaveDefaults.mat'')')


