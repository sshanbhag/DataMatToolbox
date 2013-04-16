function varargout = shufflecell(rawcell)
%------------------------------------------------------------------------
% shuffles cell across both dimensions
%	rawcell is an {N, M} cell array 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%
%------------------------------------------------------------------------
% Created: 2 April, 2013 (SJS)
%
% Revisions:
%------------------------------------------------------------------------

% get rows, cols of rawcell
[N, M] = size(rawcell);

% reshape the cell array into {N*M, 1} column vector
colcell = reshape(rawcell, N*M, 1);

% compute random permutation to shuffle indices
shuffleindices = randperm(N*M);

% shuffle and reshape to original size
varargout{1} = reshape(colcell(shuffleindices), N, M);

if nargout == 2
	varargout{2} = shuffleindices;
end