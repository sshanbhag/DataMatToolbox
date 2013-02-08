classdef Tsub < T
	properties
		cstr
		dnum
		unusedsub1
		unusedsub2
	end
	methods	
		function obj = Tsub(varargin)
			obj.cstr = 'Tsubobj';
			obj.dnum = 2001;
		end
		
		
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function varargout = getmatchproperties(obj)			
			%-----------------------------------------------
			% list the properties to match here.  
			% property_type char vector has 'n' for numeric type, 'c' for char
			%-----------------------------------------------
			property_names = {	'cstr', 'dnum'	};
			property_types = { 'char', 'num' };
			% get the superclass names and types
			[super_names, super_types] = getmatchproperties@T(obj);
			% append to the ones for this class
			property_names = [super_names property_names];
			property_types = [super_types property_types];
			
			% get property values
			if any( nargout == (0:3))
				varargout{1} = property_names;
			end
			% get property names for 2nd output arg
			if any(nargout == (2:3))
				varargout{2} = property_types;
			end
			% get property type for 3rd output arg
			if nargout == 3
				varargout{3} = cell(length(property_names), 1);
				for n = 1:length(property_names)
					varargout{1}{n} = obj.(property_names{n});
				end
			end
		end
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
				
		
		
		
		
		
		
		
	end	% End of methods
end	% End of classdef
