classdef T < handle
	properties
		a
		b
	end
	
	methods
		
		function obj = T(varargin)
			if nargin == 0
				return
			end
		end
		
		function [mval, mcomp] = match(obj, obj2)
			% first, check that classes are the same
			if ~strcmpi(class(obj), class(obj2))
				mval = 0;
				mcomp = 0;
			else
				np = 2;
				mcomp = zeros(np, 1);
				if obj.a == obj2.a
					mcomp(1) = 1;
				end
				if obj.b == obj2.b
					mcomp(2) = 1;
				end
				if all(mcomp)
					mval = 1;
				else
					mval = 0;
				end
			end
			
		end
		
	end	% End of methods
end	% End of classdef

