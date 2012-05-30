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
% 			string						raw string from text file (¿redundant?)
% 			M
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
		string
		StimulusTypeR
		StimulusTypeL
		wavFilesR
		wavFilesL
		
	% end of properties
	end
	
	%%
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods	
		
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
				obj.parseString;
			end
		end		%Marker
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
				obj.(MARKER_TAGS{m})
			end

		end	%parseString
		%---------------------------------------------------------------------	
		%---------------------------------------------------------------------	
		
	% End of methods
	end
end


