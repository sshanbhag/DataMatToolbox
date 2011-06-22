function [outSingle, outMatrix] = cellcmp(A, B)
%------------------------------------------------------------------------
%[outSingle, outMatrix] = cellcmp(A, B)
%------------------------------------------------------------------------
% given two equivalently-size cell arrays, A, and B, cellcmp will 
% compare each element and return the overall comparison result
% (1 if equal, 0 if not) in outSingle.  
%
% outMatrix will the element-by-element comparison results
%
% Limitations:
% 	- doesn't work if there are structures in the the cells
% 	- doesn't work if there are cells within cells
% 	- doesn't work if there are non-numeric arrays within cells
% 	
% 	but it WILL work if there are numeric arrays as cell elements.  note,
% 	however, that only  a single value for comparing the array elements
% 	will be returned
%
%------------------------------------------------------------------------
% Input Arguments:
% 		A				cell array of size {N, M}
% 		B				cell array of size {N, M}
%
% Output Arguments:
%		outSingle	1 if all elements of A are identical to B,
% 						0 if 1 or more elements are different
% 		outMatrix	element-by-element comparison matrix, size [N, M]
% 
%------------------------------------------------------------------------
% See also: strcmp
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 21 June, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: - recurse to handle cells?  structs? char arrays?
%------------------------------------------------------------------------

% first make sure input cells are same size
sizeA = size(A);
sizeB = size(B);

% if not same size, throw an error
if ~all(sizeA == sizeB)
	error('%s: cell A and cell B must be same size', mfilename);
end


outMatrix = zeros(sizeA);


% now, go element by element and compare values

% loop through rows
for n = 1:sizeA(1)
	% loop through cols
	for m = 1:sizeA(2)
		
		% need to check to see if this element is a string
		ischarA = ischar(A{n, m});
		ischarB = ischar(B{n, m});
		
		if (ischarA && ischarB)
			% both elements are strings, so use strcmp to compare
			outMatrix(n, m) = strcmp(A{n, m}, B{n, m});
		
		elseif (~ischarA && ~ischarB)
			% both elements are not strings, so compare with ==
			% (we use global_compare() in order to handle vectors properly)
			outMatrix(n, m) = global_compare(A{n, m}, B{n, m});
			
		else
			% these two elements cannot be equal, since
			% one is a string and the other is not - which one is not
			% important for this operation
			outMatrix(n, m) = 0;
		end
	end
end

outSingle = all(all(outMatrix));



function out = global_compare(a, b);
	tmp = (a == b);
	while ~all(size(tmp) == [1 1])
		tmp = all(tmp);
	end
	out = tmp;






