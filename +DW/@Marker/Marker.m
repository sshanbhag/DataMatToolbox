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

classdef (ConstructOnLoad = true) Marker
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
		%	Marker(<fileName>) 
		%	Constructor method
		%---------------------------------------------------------------------	

			%parse input and verify
			obj.string = '';
			if nargin == 1
				obj.string = varargin{1};
			end
		end		%Marker
		%---------------------------------------------------------------------	
		%---------------------------------------------------------------------	

		%---------------------------------------------------------------------	
		%---------------------------------------------------------------------	
		function parseString(obj)
			% loop through markers (in M() struct array), pulling out text and value
			for m = 1:MARKER_NMARKERS
				if any(strcmp(MARKER_TYPES{m}, {'int', 'float', 'double'}))
					obj.(MARKER_TAGS{m}) = str2num(string{m});
				elseif strcmp(MARKER_TYPES{m}, 'char')
					obj.(MARKER_TAGS{m}) = string{m};
				else
					error('%s: undefined marker type %s', mfilename, MARKER_TYPES{n});
				end
			end

		end	%parseString
		%---------------------------------------------------------------------	
		%---------------------------------------------------------------------	
		
	% End of methods
	end
end


