classdef  (ConstructOnLoad = true) test < handle
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	properties (SetAccess = public, GetAccess = public)
		file
		path
	end
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function obj = test(varargin)
			obj.file = '';
			obj.path = '';
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.file(obj, val)
			length(val)
			% check if no filename or path provided
			if isempty(val)
				obj.file = '';
			else
				obj.file = val;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function set.path(obj, val)
			% check if no filename or path provided
			if isempty(val)			
				obj.path = '';
			else
				% user provided path string
				obj.path = val;
			end
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

	end		% END METHODS
end		% END test CLASSDEF
