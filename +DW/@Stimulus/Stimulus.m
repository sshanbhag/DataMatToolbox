%-----------------------------------------------------------------------------
% Stimulus.m
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW package
% Class Definition
%-----------------------------------------------------------------------------
%	Properties:
% 		Channel						'R', 'L', 'B'
%		Amplitude					output signal amplitude (mV)
% 		Attenuation
% 		TimeShift
% 		RampUp
% 		HoldTime
% 		RampDown
% 		OutputTimestamp
% 		OutputTimeWithDelay
% 		FixedDelay
%-----------------------------------------------------------------------------
% See also: Wav, Tone, Noise, StimulusList, Data, Marker, Probe, Unit
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 4 June, 2012 (SJS)
%
% Revisions:
%	12 Feb 2013 (SJS): reworked, subclassed (Tone Wav Noise)
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% classdef (ConstructOnLoad = true) Stimulus < handle
classdef Stimulus < handle
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Protected Properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		Channel
		Amplitude
		Attenuation
		TimeShift
		RampUp
		HoldTime
		RampDown
		OutputTimestamp
		OutputTimeWithDelay
		FixedDelay
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
		function obj = Stimulus(varargin)
		%---------------------------------------------------------------------	
		% Stimulus
		% Constructor method
		%---------------------------------------------------------------------
		% Stimulus()	when called with no arguments, returns empty
		%					Stimulus object
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;
			%--------------------------------------------------------
			% Check inputs
			%--------------------------------------------------------
			% if no inputs, return
			if isempty(varargin)
				return
			end
			% if 1 input, need to perform a few checks
			if length(varargin) == 1
				% check if varargin{1} was directly passed in from
				% a subclass (using obj@DW.Stimulus(varargin) syntax)
				if strcmpi(inputname(1), 'varargin')
					% make local copy of varargin{1}
					tmp = varargin{1};
					if length(tmp) == 1
						% tmp is marker only, assume RIGHT channel
						obj.buildStimulusFromMarker(tmp{1}, R);
					elseif length(tmp) == 2
						% caller provided Marker and Channel as inputs
						obj.buildStimulusFromMarker(tmp{1}, tmp{2});
					else
						fprintf('%s: empty or strange inputs... ', mfilename);
					end
					clear tmp
				else
					% if no channel provided, assume channel it RIGHT
					obj.buildStimulusFromMarker(varargin{1}, R);
				end
			elseif length(varargin) == 2
				% caller provided Marker and Channel as inputs
				obj.buildStimulusFromMarker(varargin{1}, varargin{2});
			else
				error('%s: invalid input args', mfilename);
			end

		end		% END Stimulus constructor
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		% General Methods
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function obj = buildStimulusFromMarker(obj, M, C)
		%---------------------------------------------------------------------
		% Stimulus.buildStimulusFromMarker(M, C)
		%---------------------------------------------------------------------
		% given marker object M and channel C (L == 1, R == 2), assign values
		% from marker to appropriate values in Stimulus object
		%---------------------------------------------------------------------
			DW.DataWaveDefaults;
			markerbase = {	'Amplitude', 'Attenuation', 'TimeShift', ...
							'RampUp', 'HoldTime', 'RampDown', ...
							'FixedDelay', 'OutputTimestamp', 'OutputTimeWithDelay' };
			% set channel, and marker tag character (to postpend to
			% markerbase string to properly identify Marker property name)
			if C == L
				tagchar = 'L';
				obj.Channel = L;
			else
				tagchar = 'R';
				obj.Channel = R;
			end
			% assign values to tags in Stimulus object
			for n = 1:length(markerbase)
				obj.(markerbase{n}) = M.([markerbase{n} tagchar]);
			end

		end	% END buildStimulusFromMarker
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------

		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function varargout = getmatchproperties(obj)
		%---------------------------------------------------------------------
		% [names, type, values] = Stimulus.getmatchproperties
		%---------------------------------------------------------------------
			
			%-----------------------------------------------
			% list the properties (and type) to match here.  
			% Amplitude is left out simply because it is often used as an
			% adjustment during calibration.  
			% !!! this list may need to be adjusted in the future!!!!
			% property_type char vector has 'n' for numeric type, 'c' for char
			% this will have to be modified if properties are added 
			% subclasses will also need to act accordingly!
			%-----------------------------------------------
			property_names = {	'Channel', 'TimeShift', 'RampUp', 'HoldTime', ...
								'RampDown', 'FixedDelay'	};
			property_types = { 'num', 'num', 'num', 'num', 'num', 'num' };
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
		function [mval, mcomp] = match(obj, B)
		%---------------------------------------------------------------------
		% [mval, mcomp] = Stimulus.match(stimobj)
		%---------------------------------------------------------------------
		
			%-----------------------------------------------------------
			% get the property names and types to match 
			%-----------------------------------------------------------
			[mProp, mType] = obj.getmatchproperties;
			% get lengths
			np = length(mProp);	
			nB = length(B);
			% preallocate things
			mval = zeros(nB, 1);
			mcomp = cell(nB, 1);

			%-----------------------------------------------------------
			% loop through the input B (can be a singleton or a vector!)
			%-----------------------------------------------------------
			for b = 1:nB;
				% local copy of element in B (in order to allow for
				% arrays or cells of data)
				if iscell(B)
					Btmp = B{b};
				else
					Btmp = B(b);
				end
				%-----------------------------------------------
				% first, check that classes are the same
				%-----------------------------------------------
				if ~strcmpi(class(obj), class(Btmp))
					% if they're not the same, assign zeros
					mval(b) = 0;
					mcomp{b} = 0;
				else
					%-----------------------------------------------
					% since they are the same, compare values
					%-----------------------------------------------
					% preallocate zeros vector for item comparisons
					mcomp{b} = zeros(np, 1);
					% loop through properties to compare
					for n = 1:np
						% compare using appropriate method for property type
						if strcmpi(mType{n}, 'num')
							% number!
							mcomp{b}(n) = ( obj.(mProp{n}) == Btmp.(mProp{n}) );
						elseif strcmpi(mType{n}, 'char')
							% char (or string)!
							mcomp{b}(n) = strcmp(obj.(mProp{n}), Btmp.(mProp{n}));
						else
							% WTF????
							error('%s: unknown type %s!!!!', mfilename, mType{n});
						end
					end
					if all(mcomp{b})
						% if they're all 1, they match
						mval(b) = 1;
					else
						% if not, they don't and mval(b) is 0
						mval(b) = 0;
					end
				end	% END if
			end	% END b
		end	% END match FUNCTION
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
end	% End of classdef






