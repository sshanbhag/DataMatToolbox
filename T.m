classdef T < handle
	properties
		anum
		bstr
		unused1
		unused2
	end
	
	methods
		
		function obj = T(varargin)
			obj.anum = 1001;
			obj.bstr = 'Tobj';
		end
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function varargout = getmatchproperties(obj)			
			%-----------------------------------------------
			% list the properties to match here.  
			% property_type char vector has 'n' for numeric type, 'c' for char
			%-----------------------------------------------
			property_names = {	'anum', 'bstr'	};
			property_types = { 'num',  'char' };
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
		
		%---------------------------------------------------------------------
		%---------------------------------------------------------------------
		function [mval, mcomp] = match(obj, B)
		%---------------------------------------------------------------------
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
				%-----------------------------------------------
				% first, check that classes are the same
				%-----------------------------------------------
				if ~strcmpi(class(obj), class(B(b)))
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
							mcomp{b}(n) = ( obj.(mProp{n}) == B(b).(mProp{n}) );
						elseif strcmpi(mType{n}, 'char')
							% char (or string)!
							mcomp{b}(n) = strcmp(obj.(mProp{n}), B(b).(mProp{n}));
						else
							% WTF????
							error('%s: unknown type %s!!!!', mfilename, mType{b});
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
		
	end	% End of methods
end	% End of classdef

