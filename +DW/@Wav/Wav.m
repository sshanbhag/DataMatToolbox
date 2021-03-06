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
%	? day, 2013 (SJS): changed to use DW.fileparts so that PC behaviour for
%		fileparts is used (path works successfully on all platforms)
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
			DW.DataWaveDefaults;
			%--------------------------------------------------------
			% call superclass constructor
			%--------------------------------------------------------
			obj = obj@DW.Stimulus(varargin);
			%--------------------------------------------------------
			% Noise-specific init
			%--------------------------------------------------------
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
			DW.DataWaveDefaults;
			% set frequency
			if Channel == R
				[obj.Filepath, tmpn, tmpx] = ...
												DW.fileparts(Marker.WavFilenameR, 'PC');
			else
				[obj.Filepath, tmpn, tmpx] = ...
												DW.fileparts(Marker.WavFilenameL, 'PC');
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
		%---------------------------------------------------------------------
		function varargout = getmatchproperties(obj)
		%---------------------------------------------------------------------
		% [names, types, values] = Wav.getmatchproperties
		%---------------------------------------------------------------------
		% Method to return properties used by Stimulus.match method
		% (superclass).  
		% 	names is N X 1 cell vector of property names
		% 	types is N X 1 cell vector of property types
		% 						'num' indicates numeric type
		% 						'char' indicates character or string type
		% 	values is N X 1 cell vector of valus of properties listed in names
		%---------------------------------------------------------------------
			
			%-----------------------------------------------
			% list the properties to match here.  
			%-----------------------------------------------
			property_names = {	'Filename'	};
			property_types = { 'char' };
			nprop = length(property_names);
			%-----------------------------------------------
			% get property values
			%-----------------------------------------------
			if any( nargout == (0:3))
				varargout{1} = property_names;
			end
			%-----------------------------------------------
			% get property names for 2nd output arg
			%-----------------------------------------------
			if any(nargout == (2:3))
				varargout{2} = property_types;
			end
			%-----------------------------------------------
			% get property type for 3rd output arg
			%-----------------------------------------------
			if nargout == 3
				varargout{3} = cell(length(property_names), 1);
				for n = 1:length(property_names)
					varargout{1}{n} = obj.(property_names{n});
				end
			end
		end	% END getmatchproperties
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






