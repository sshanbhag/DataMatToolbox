function [uniqueVals, uniqueIndices, NuniqueVals] = findUniqueValues(values_to_search)
%-----------------------------------------------------------------------------
%[uniqueVals, uniqueIndices, NuniqueVals] = findUniqueValues(values_to_search)
%-----------------------------------------------------------------------------
% Finds unique values in list of values_to_search
% 						
%-----------------------------------------------------------------------------
% Input Arguments:
%	values_to_search
%
% Output Arguments:
%	uniqueVals		unique values of values_to_search array
%	uniqueIndices	indices where each of the unique values are found in values_to_search
%	NuniqueVals		# of unique values in values_to_search
%-----------------------------------------------------------------------------
% See Also: findUniqueStrings, findUniqueCellRows, unique
%-----------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 3 July, 2011 (SJS)
% 	- modified from findUniqueText.m
%
% Revisions:
%	6 July, 2011 (SJS): updated documentation
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% Initialize some arrays and variables
Nvalues = length(values_to_search);

testIndex = 1;
runFlag = 1;


% find the first instance of each unique value
uniqueVals = unique(values_to_search, 'first');

NuniqueVals = length(uniqueVals);

if ~NuniqueVals
	uniqueIndices = {};	
	return
else
	uniqueIndices = cell(NuniqueVals, 1);
end

for n = 1:NuniqueVals
	uniqueIndices{n} = find(uniqueVals(n) == values_to_search);
end
