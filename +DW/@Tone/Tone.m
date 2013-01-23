%-----------------------------------------------------------------------------
% Tone Class
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW package
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
%-----------------------------------------------------------------------------
% See also: Stimulus, Noise, Wav, DWdata, Marker, Probe, Unit
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 18 January, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*** important: superclass declaration must have package name (DW) prepended!
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
classdef (ConstructOnLoad = true) Tone < DW.Stimulus
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	properties	(SetAccess = protected)
		Freq
		Phase
	end	% end of protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	methods	
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Constructor Method
		%---------------------------------------------------------------------
		function obj = Tone(varargin)
		%---------------------------------------------------------------------	
		% Tone < DW.Stimulus
		% Constructor method
		%---------------------------------------------------------------------
		% Tone()	when called with no arguments, returns empty
		%			Tone object
		%---------------------------------------------------------------------
			DataWaveDefaults;
			%--------------------------------------------------------
			% call superclass constructor
			%--------------------------------------------------------
			obj = obj@DW.Stimulus(varargin);
			%--------------------------------------------------------
			% Noise-specific init
			%--------------------------------------------------------
			% set Type (def. in Stimulus) to NOISE
			obj.Type = 'TONE';
			if isempty(varargin)
				return
			elseif length(varargin) == 1
				% if no c
				C = R;
			else
				C = varargin{2};
			end
			obj.setValsFromMarker(varargin{1}, C);
		end	% END Tone constructor
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% General Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function obj = setValsFromMarker(obj, Marker, Channel)
		%---------------------------------------------------------------------
			DataWaveDefaults;
			% set frequency
			if Channel == R
				obj.Freq = Marker.ToneFreqR;
				obj.Phase = Marker.PhaseDegR;
			else
				obj.Freq = Marker.ToneFreqL;
				obj.Phase = Marker.PhaseDegL;
			end
		end	% END setValsFromMarker
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% Overloaded Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% set/get Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
	
	end	% End of methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
end	% End of classdef
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************




