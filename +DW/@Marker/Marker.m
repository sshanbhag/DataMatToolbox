%-----------------------------------------------------------------------------
% Marker.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% Class Definition
%-----------------------------------------------------------------------------
% 			string						raw string from text file
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
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		string
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
		M
		StimulusTypeR
		StimulusTypeL
		wavFilesR
		wavFilesL
	end
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods
	
		
	%% ------------------------------------------------------------------------
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
	
%% ------------------------------------------------------------------------
	function parseString(obj)
	
	end	%parseString
		
	% End of methods
	end
end


