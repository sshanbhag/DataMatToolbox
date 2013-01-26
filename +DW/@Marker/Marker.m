%-----------------------------------------------------------------------------
% Marker.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% Class Definition
%-----------------------------------------------------------------------------
% 	marker "tags" from DataWave/Neuroshare events :
% 		SoundTypeR
% 		AttenuationR
% 		WavFilenameR
% 		BBNlowerFreqR
% 		BBNupperFreqR
% 		AmplitudeR
% 		TimeShiftR
% 		RampUpR
% 		HoldTimeR
% 		RampDownR
% 		OutputTimestampR
% 		OutputTimeWithDelayR
% 		FixedDelayR
% 		PA5idR
% 		ToneFreqR
% 		PhaseDegR
% 		OutputFileR
% 		SoundTypeL
% 		AttenuationL
% 		WavFilenameL
% 		BBNlowerFreqL
% 		BBNupperFreqL
% 		AmplitudeL
% 		TimeShiftL
% 		RampUpL
% 		HoldTimeL
% 		RampDownL
% 		OutputTimestampL
% 		OutputTimeWithDelayL
% 		FixedDelayL
% 		PA5idL
% 		ToneFreqL
% 		PhaseDegL
% 		OutputFileL
% 
%	Additional Parameters:
% 		Timestamp
% 		ID
% 		string
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
%	14 Jan 2013 (SJS)
%	 -	modified to use DDF event data from NeuroShare
%-----------------------------------------------------------------------------
% TO DO:
%
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
classdef (ConstructOnLoad = true) Marker < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		% tags
		SoundTypeR
		AttenuationR
		WavFilenameR
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
		ToneFreqR
		PhaseDegR
		OutputFileR
		SoundTypeL
		AttenuationL
		WavFilenameL
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
		ToneFreqL
		PhaseDegL
		OutputFileL
		Timestamp
		ID
		string
	end	% end of protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Constructor Method
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function obj = Marker(varargin)
		%---------------------------------------------------------------------	
		% Marker
		%---------------------------------------------------------------------	
		% Constructor method
		%---------------------------------------------------------------------
		% Marker(dwstring)	where dwstring is a row of text, in cell form, 
		%							from a datawave text file, will create an
		%							initialized/parsed Marker object
		%							
		%							when called with no arguments, returns empty
		%							Marker object
		%---------------------------------------------------------------------
			nargs = length(varargin);	% check input args
			% if no args provided, simply return
			if nargs == 0
				return
			end
			% mimick old behavior for strings - parse input and verify
			% !!! this should probably disappear (SJS 17Jan2013)
			if nargs == 1
				obj.string = '';
				obj.string = varargin{1};
				obj.initMarkerFromString;
				return
			end
			% otherwise, parse inputs args
			argc = 1;
			while argc <= nargs
				switch upper(varargin{argc})
					case 'EVENT_STRING'
						obj.setValuesFromEventList(varargin{argc+1});
						argc = argc + 2;
					otherwise
						error('%s: unknown option %s', mfilename, varargin{argc})
				end
			end
		end	% END Marker
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Comparison Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function [out, outcmp, outmat] = sameStimulus(obj, objA, channel)
		%---------------------------------------------------------------------
		% [out, outcmp, outmat] = Marker.sameStimulus(objA, channel)
		%---------------------------------------------------------------------
		% compares the stimulus-specific parameters from Marker to those of
		% another vector, objA.  
		% These parameters are:
		% 	SoundType
		% 	Attenuation
		% 	WavFilename
		% 	BBNlowerFreq
		% 	BBNupperFreq
		% 	Amplitude
		% 	TimeShift
		% 	RampUp
		% 	HoldTime
		% 	RampDown
		% 	FixedDelay
		% 	ToneFreq
		% 	PhaseDeg
		%
		% if channel == 'L', the values will be from the Left channel (e.g.,
		% SoundTypeL, AttenuationL, etc.)
		% channel == 'R' will compare Right channel values.
		% if channel is not provided, or if it equals 'B' (both), the 
		% R and L channel values will be compared
		%
		% out == 1 indicates stimuli are the same
		% out == 0 indicates different stimuli
		% outcmp is a list of comparisons between the values
		% outmat is a {n, 2} list of the values that were compared
		%---------------------------------------------------------------------
			DataWaveDefaults;	% load defaults
			[out, outcmp] = cellcmp(obj.getStimulus(channel), ...
											objA.getStimulus(channel));
			if nargout == 3
				outmat = [obj.getStimulus(channel) objA.getStimulus(channel)];
			end
		end	% END Marker/eq
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% General Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function out = getStimulus(obj, channel)
		%---------------------------------------------------------------------
		% gets the stimulus-specific parameters from Marker, returns
		% a cell vector with list of values:
		% 	SoundType
		% 	Attenuation
		% 	WavFilename
		% 	BBNlowerFreq
		% 	BBNupperFreq
		% 	Amplitude
		% 	TimeShift
		% 	RampUp
		% 	HoldTime
		% 	RampDown
		% 	FixedDelay
		% 	ToneFreq
		% 	PhaseDeg
		%
		% if channel == 'L', the values will be from the Left channel (e.g.,
		% SoundTypeL, AttenuationL, etc.)
		% channel == 'R' will return Right channel values.
		% if channel is not provided, or if it equals 'B' (both), the 
		% output vector will have the R channel values followed by the L 
		% channel values
		%---------------------------------------------------------------------
			DataWaveDefaults;	% load defaults
			% use different indices into tags depending on channel
			if nargin == 1
				% get all tags if channel not provided
				tags = MARKER_TAGS([RCOMP_INDEX LCOMP_INDEX]);
			else
				if channel == 'L'
					tags = MARKER_TAGS(LCOMP_INDEX);
				elseif channel == 'R' 
					tags = MARKER_TAGS(RCOMP_INDEX);
				else
					tags = MARKER_TAGS([RCOMP_INDEX LCOMP_INDEX]);
				end
			end
			% get # of tags, preallocate output array
			nout = length(tags);
			out = cell(nout, 1);
			% get values
			for n = 1:nout
				out{n} = obj.(tags{n});
			end
		end	% END getStimulus
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function setValuesFromEventList(obj, Elist)
		%---------------------------------------------------------------------
		% Marker.setValuesFromEventList(Elist)
		%---------------------------------------------------------------------
		% it is assumed that Elist is a cell array of values (pre-converted
		% from strings) for the Marker parameters (listed in MARKER_TAGS)
		%---------------------------------------------------------------------
			DataWaveDefaults;		% load defaults
			% check that Elist is a cell
			if ~iscell(Elist)
				error('%s: list of values must be in cell format', mfilename);
			end
			if length(Elist) ~= MARKER_NMARKERS
				error('%s: mismatch in # of fields in Elist', mfilename);
			end
			% assign values, using list of MARKER_TAGS to index into
			% the Marker parameters
			for n = 1:MARKER_NMARKERS
				obj.(MARKER_TAGS{n}) = Elist{n}; %#ok<*USENS>
			end
		end
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
				
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
			DataWaveDefaults;	% load defaults
			if length(varargin) == 1
				obj.string = varargin{1};
			end
			% loop through markers, pulling out text and value
			for m = 1:MARKER_NMARKERS
				% check if current marker is a number
				if any(strcmp(MARKER_TYPES{m}, {'int', 'float', 'double'}))
					% if number, store value in vector 
					obj.(MARKER_TAGS{m}) = str2num(obj.string{m}); %#ok<ST2NM>
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
		%---------------------------------------------------------------------
			DataWaveDefaults;	% load defaults
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
			DataWaveDefaults;				% load defaults
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

		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Set/Get methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function setTimestamp(obj, val)
			if length(val) == 1
				obj.Timestamp = val;
				return
			else
				if all(val > 0)
					% check if they are equal
					if val(1) == val(2)
						% if so, use the right one
						obj.Timestamp = val(1);
					else
						% otherwise, store the first one (in time)
						obj.Timestamp = min(val);
					end
				else
					% store timestamp that is > 0
					obj.Timestamp = val(val>0);
				end
			end			
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function setID(obj, val)
			if ~isempty(val)
				obj.ID = val;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

%{
		%------------------------------------------------------------------------
		% set N (number?)
		%------------------------------------------------------------------------
		function set.N(obj, val)
			if isempty(val)
				% check if no value - if so, set to empty
				obj.N = [];
				return
			elseif isnumeric(val)
				obj.N = val;
				return
			end
		end
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.Timestamp(obj, val)
			if isempty(val)
				return;
			elseif isnumeric(val)
				obj.Timestamp = val;
			end
			return
		end
		%------------------------------------------------------------------------
%}
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Overloaded Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function out = eq(objA, objB)
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
			DataWaveDefaults;	% load defaults
			% find fields to compare
			numfields = find( strcmp('int', MARKER_TYPES) | ...
									strcmp('float', MARKER_TYPES) );
			charfields = find( strcmp('char', MARKER_TYPES) );
			% preallocate comparison arrays
			numcomp = zeros(1, length(numfields));
			charcomp = zeros(1, length(charfields));
			% compare numeric values
			for n = 1:length(numfields)
				t = numfields(n);
				numcomp(n) = all(objA.(MARKER_TAGS{t}) == objB.(MARKER_TAGS{t}));
			end
			% compare character values
			for c = 1:length(charfields)
				t = charfields(c);
				charcomp(c) = strcmp( objA.(MARKER_TAGS{t}), objB.(MARKER_TAGS{t}) );
			end
			% all check?
			out = all([numcomp charcomp]);
		end	% END Marker/eq
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
	end	% End of methods
end


