%-----------------------------------------------------------------------------
% Noise Class
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW package
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
%-----------------------------------------------------------------------------
% See also: Stimulus, Tone, Wav, DWdata, Marker, Probe, Unit
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
classdef (ConstructOnLoad = true) Noise < DW.Stimulus
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Protected Properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		LowerFreq
		UpperFreq
	end	% end of properties
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
		function obj = Noise(varargin)
		%---------------------------------------------------------------------	
		% Noise < DW.Stimulus
		% Constructor method
		%---------------------------------------------------------------------
		% Noise()	when called with no arguments, returns empty
		%				Noise object
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
			obj.Type = 'NOISE';
			if isempty(varargin)
				return
			elseif length(varargin) == 1
				% if no c
				C = R;
			else
				C = varargin{2};
			end
			obj.setValsFromMarker(varargin{1}, C);
		end	% END Noise constructor
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
				obj.LowerFreq = Marker.BBNlowerFreqR;
				obj.UpperFreq = Marker.BBNupperFreqR;
			else
				obj.LowerFreq = Marker.BBNlowerFreqL;
				obj.UpperFreq = Marker.BBNupperFreqL;
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

