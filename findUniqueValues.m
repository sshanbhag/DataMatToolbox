function [uniqueVals, uniqueIndices, NuniqueVals] = findUniqueValues(values_to_search)
%-----------------------------------------------------------------------------
%[uniqueVals, uniqueIndices, NuniqueVals] = findUniqueText(values_to_search)
%-----------------------------------------------------------------------------
%	algorithm: 
% 		1) obtain target text string (searchstr) from list values_to_search
% 		2) find locations of searchstr in full list of values_to_search
% 		3) see if this string already exists in the uniqueVals list
% 				If yes:
%					- increment the testIndex so that the next value in
% 					  testvals becomes the new searchstr
% 				If no:
% 					- increment NuniqueVals count
% 					- store the current searchstr in uniqueVals
% 					- store indices of searchstr in values_to_search in the uniqueIndices list
% 					- rebuild the testvals list by removing instances of the current
% 					  searchstr
% 						
%-----------------------------------------------------------------------------
% NuniqueVals		# of unique strings in values_to_search
% uniqueVals		unique values of strings in values_to_search, cell array
% uniqueIndices	indices where each of the unique strings are found in values_to_search
%-----------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 3 July, 2011 (SJS)
% 	- modified from findUniqueText.m
%
% Revisions:
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
	uniqueIndices{n} = find(uniqueVals(n) == values_to_search)
end
