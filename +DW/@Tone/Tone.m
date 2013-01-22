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
	% Define properties
	%------------------------------------------------------------------------
	properties
		Freq
		Phase
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
		function obj = Tone(varargin)
		%---------------------------------------------------------------------	
		% Tone < DW.Stimulus
		% Constructor method
		%---------------------------------------------------------------------
		% Tone()	when called with no arguments, returns empty
		%			Tone object
		%---------------------------------------------------------------------

			%--------------------------------------------------------
			%parse input and verify
			%--------------------------------------------------------
			obj.Type = 'TONE';
			obj.Freq = [];
			obj.Phase = [];
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
		function obj = buildStimulusFromMarkers(obj)
		%---------------------------------------------------------------------	


		end	% END buildStimulusFromMarkers
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




