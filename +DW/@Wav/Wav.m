%-----------------------------------------------------------------------------
% Wav Class
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW package
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
%-----------------------------------------------------------------------------
% See also: Stimulus, Noise, Tone, DWdata, Marker, Probe, Unit
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
% classdef (ConstructOnLoad = true) Wav < DW.Stimulus
classdef Wav < DW.Stimulus
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Protected Properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		Filename
		Filepath
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
		function obj = Wav(varargin)
		%---------------------------------------------------------------------	
		% Wav < DW.Stimulus
		% Constructor method
		%---------------------------------------------------------------------
		% Wav()	when called with no arguments, returns empty
		%				Wav object
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
			obj.Type = 'WAV';
			if isempty(varargin)
				return
			elseif length(varargin) == 1
				% if no c
				C = R;
			else
				C = varargin{2};
			end
			obj.setValsFromMarker(varargin{1}, C);
		end	% END Wav constructor
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
				[obj.Filepath, tmpn, tmpx] = fileparts(Marker.WavFilenameR);
			else
				[obj.Filepath, tmpn, tmpx] = fileparts(Marker.WavFilenameL);
			end
			obj.Filename = [tmpn tmpx];
			clear tmpn tmpx
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






