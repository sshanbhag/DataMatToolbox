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
% classdef (ConstructOnLoad = true) Noise < DW.Stimulus
classdef Noise < DW.Stimulus
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
		function varargout = getmatchproperties(obj)
		%---------------------------------------------------------------------
		% [names, types, values] = Noise.getmatchproperties
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
			property_names = {	'LowerFreq', 'UpperFreq'	};
			property_types = {	'num', 'num'	};		% c = char, n = num
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

