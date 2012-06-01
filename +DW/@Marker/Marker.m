%-----------------------------------------------------------------------------
% Marker.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% Class Definition
%-----------------------------------------------------------------------------
%		marker "tags":
% 			Timestamp					time of event, in microseconds
% 			id								?
% 			OutputFilename				output data file name
% 			AttenuationR				R channel attenuation setting (dB)
% 			BBNlowerFreqR				if noise stimulus, 
% 			BBNupperFreqR
% 			AmplitudeR
% 			TimeShiftR
% 			RampUpR
% 			HoldTimeR
% 			RampDownR
% 			OutputTimestampR
% 			OutputTimeWithDelayR
% 			FixedDelayR
% 			PA5idR
% 			WavFilenameR
% 			ToneFreqR
% 			PhaseDegR
% 			AttenuationL
% 			BBNlowerFreqL
% 			BBNupperFreqL
% 			AmplitudeL
% 			TimeShiftL
% 			RampUpL
% 			HoldTimeL
% 			RampDownL
% 			OutputTimestampL
% 			OutputTimeWithDelayL
% 			FixedDelayL
% 			PA5idL
% 			WavFilenameL
% 			ToneFreqL
% 			PhaseDegL
%			
%		Added values
%			N
% 			string						raw string from text file (¿redundant?)
% 			StimulusTypeR
% 			StimulusTypeL
% 			wavFilesR
% 			wavFilesL
%-----------------------------------------------------------------------------
% See also: DWdata, loadDWfile (function)
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 17 May, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%
%-----------------------------------------------------------------------------

classdef (ConstructOnLoad = true) Marker < handle
	%%
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		% tags
		Timestamp
		id
		OutputFilename
		AttenuationR
		BBNlowerFreqR
		BBNupperFreqR
		AmplitudeR
		TimeShiftR
		RampUpR
		HoldTimeR
		RampDownR
		OutputTimestampR
		OutputTimeWithDelayR
		FixedDelayR
		PA5idR
		WavFilenameR
		ToneFreqR
		PhaseDegR
		AttenuationL
		BBNlowerFreqL
		BBNupperFreqL
		AmplitudeL
		TimeShiftL
		RampUpL
		HoldTimeL
		RampDownL
		OutputTimestampL
		OutputTimeWithDelayL
		FixedDelayL
		PA5idL
		WavFilenameL
		ToneFreqL
		PhaseDegL
		% other
		StimulusTypeR
		StimulusTypeL
		
	% end of properties
	end

	properties
		N
		string
	end

	
	%%
	%%
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods	
		
		%% Constructor Method
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function obj = Marker(varargin)
		%---------------------------------------------------------------------	
		% Marker
		% Constructor method
		%---------------------------------------------------------------------
		% Marker(dwstring)	where dwstring is a row of text, in cell form, 
		%							from a datawave text file, will create an
		%							initialized/parsed Marker object
		%							
		%							when called with no arguments, returns empty
		%							Marker object
		%---------------------------------------------------------------------

			%--------------------------------------------------------
			%parse input and verify
			%--------------------------------------------------------
			obj.string = '';
			
			if nargin == 1
				obj.string = varargin{1};
				obj.initMarkerFromString;
			end
		end		%Marker
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%% Overloaded Methods

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function out = eq(objA, objB)
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
			% define fields to compare
			numfields = {	...
								'AttenuationR', ...
								'BBNlowerFreqR', ...
								'BBNupperFreqR', ...
								'AmplitudeR', ...
								'TimeShiftR', ...
								'RampUpR', ...
								'HoldTimeR', ...
								'RampDownR', ...
								'FixedDelayR', ...
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
								'FixedDelayL' ...
							};
			charfields = {	...
								'WavFilenameR', ...
								'WavFilenameL', ...
								'ToneFreqL', ...
								'PhaseDegL' ...
							};

			numcomp = zeros(size(numfields));
			charcomp = zeros(size(charfields));
			
			for n = 1:length(numfields)
				numcomp(n) = all(objA.(numfields{n}) == objB.(numfields{n}));
			end
			
			for c = 1:length(charfields)
				charcomp(c) = strcmp( objA.(charfields{c}), objB.(charfields{c}) );
			end
			
			out = all([numcomp charcomp]);

		end	% END Marker/eq
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		
		%% General Methods
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function initMarkerFromString(obj, varargin)
		%---------------------------------------------------------------------
		% initMarkerFromString(obj)
		%---------------------------------------------------------------------
		% 
		%---------------------------------------------------------------------
			if length(varargin) == 1
				obj.string = varargin{1};
			end
			if ~isempty(obj.string)
				obj.parseString;
				obj.identifyStimulus;
			else
				warning('empty string')
			end
		end
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function parseString(obj, varargin)
		%---------------------------------------------------------------------
		% parseString(obj)
		%---------------------------------------------------------------------
		% parses string into individual properties for Marker object
		%---------------------------------------------------------------------

			% load defaults
			DataWaveDefaults;
			
			if length(varargin) == 1
				obj.string = varargin{1};
			end
			
			% loop through markers, pulling out text and value
			for m = 1:MARKER_NMARKERS
				% check if current marker is a number
				if any(strcmp(MARKER_TYPES{m}, {'int', 'float', 'double'}))
					% if number, store value in vector 
					obj.(MARKER_TAGS{m}) = str2num(obj.string{m});
				elseif strcmp(MARKER_TYPES{m}, 'char')
					% otherwise, if string, store in 1-D cell array
					obj.(MARKER_TAGS{m}) = obj.string{m};
				else
					% unknown type, throw error
					error('%s: undefined marker type %s', mfilename, MARKER_TYPES{n});
				end
			end

		end	%parseString
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function identifyStimulus(obj, varargin)
		%---------------------------------------------------------------------
		% identifyStimulus(obj)
		%---------------------------------------------------------------------
		% 
		%---------------------------------------------------------------------

			% load defaults
			DataWaveDefaults;
			
			% Identify Stimulus types, store in obj.Marker structure
			obj.checkStimulusType('R');
			obj.checkStimulusType('L');
			
		end	%identifyStimulus
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		function StimulusType = checkStimulusType(obj, C)
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

			% build some tagnames
			wavtag = ['WavFilename' C];
			tonefreqtag = ['ToneFreq' C];
			bbnlowtag = ['BBNlowerFreq' C];
			bbnhitag = ['BBNupperFreq' C];
		
			%----------------------------------------
			% first check if wav filename exists
			%----------------------------------------
			if isempty(obj.(wavtag))
				% no wav filename, set wavFlag to 0;
				wavFlag = 0;
			else
				wavFlag = 1;
			end

			%----------------------------------------
			% check for tone
			%----------------------------------------
			if obj.(tonefreqtag) > 0
				toneFlag = 1;
			else
				toneFlag = 0;
			end

			%----------------------------------------
			% check for noise
			%----------------------------------------
			% low frequency
			if obj.(bbnlowtag) > 0
				BBNlowFlag = 1;
			else
				BBNlowFlag = 0;
			end
			% high frequency
			if obj.(bbnhitag) > 0
				BBNhiFlag = 1;
			else
				BBNhiFlag = 0;
			end

			if BBNlowFlag && BBNhiFlag
				% if both bbn low and bbn hi flags are set, broadband noise
				% stimulus was used; set noiseFlag to 1
				noiseFlag = 1;
			elseif (BBNlowFlag && ~BBNhiFlag) || (~BBNlowFlag && BBNhiFlag)
				% for some reason, only one of the noise freq. flags was set - 
				% deliver warning to user, set noiseFlag to 0;
				warning('%s: bizarre condition for stimulus - either BBNlowerF or BBNupperF is zero', ...
								mfilename);
				unknownFlag = 1;
				noiseflag = 0;
			else
				noiseFlag = 0;
			end

			%-----------------------------------------------------
			% now perform tests and determine stimulus type
			%-----------------------------------------------------
			% build vector from flags to simplify checking
			flagVector = [wavFlag toneFlag noiseFlag];

			if unknownFlag
				% stimulus type could not be determined
				StimulusType = 'UNKNOWN';
				warning('%s: StimulusType%c is unknown', mfilename, C);

			elseif ~any(flagVector)
				% no sound for this channel
				StimulusType = 'NO_SOUND';

			else
				% check if more than one flag is set
				if sum(flagVector) > 1 
					warning('%s: multiple stimuli for sweep # %d', mfilename, obj.N);
					disp(flagVector);
					multFlag = 1;
				else
					multFlag = 0;
				end

				if wavFlag
					StimulusType = 'WAVFILE';
				end
				if toneFlag
					if multFlag
						StimulusType = [StimulusType 'TONE'];
					else
						StimulusType = 'TONE';
					end
				end
				if noiseFlag
					if multFlag
						StimulusType = [StimulusType 'NOISE'];
					else
						StimulusType = 'NOISE';
					end
				end
			end
			
			%-----------------------------------------------------
			% store stimulus type in object
			%-----------------------------------------------------
			if C == 'R'
				obj.StimulusTypeR = StimulusType;
			elseif C == 'L'
				obj.StimulusTypeL = StimulusType;
			end

		end
		%------------------------------------------------------------------------

		
		%% Set/Get methods
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.N(obj, val)
			% check if no value
			if isempty(val)			
				obj.N = '';
				return
			elseif isnumeric(val)
				obj.N = val;
				return
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
	% End of methods
	end
end


