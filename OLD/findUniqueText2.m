function [uniqueText, uniqueIndices, NuniqueText] = findUniqueText2(strings_to_search)
%-----------------------------------------------------------------------------
%[uniqueText, uniqueIndices, NuniqueText] = findUniqueText(strings_to_search)
%-----------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------
% Finds unique strings within the list of strings_to_search
% 
%	algorithm: 
% 		1) obtain target text string (searchstr) from list strings_to_search
% 		2) find locations of searchstr in full list of strings_to_search
% 		3) see if this string already exists in the uniqueText list
% 				If yes:
%					- increment the testIndex so that the next value in
% 					  testvals becomes the new searchstr
% 				If no:
% 					- increment NuniqueText count
% 					- store the current searchstr in uniqueText
% 					- store indices of searchstr in strings_to_search in the 
%					  uniqueIndices list
% 					- rebuild the testvals list by removing instances of the current
% 					  searchstr
%-----------------------------------------------------------------------------
% Input Arguments:
%	strings_to_search		cell array of strings
%
% Output Arguments:
% 	uniqueText		unique values of strings in strings_to_search, cell array
% 	uniqueIndices	indices where each of the unique strings are found 
%						in strings_to_search
% 	NuniqueText		# of unique strings in strings_to_search
%-----------------------------------------------------------------------------
% See Also: findUniqueStrings, findUniqueCellRows, unique
%-----------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 13 June, 2011 (SJS)
% 	- uses code snipped from loadDWStimData.m
%
% Revisions:
%	6 July, 2011 (SJS): updated documentation
%	18 Jan 2013 (SJS): minor tweaks, updated some comments/docs
%	18 Jan 2012 (SJS): 
%------------------------------------------------------------------------
% TO DO: any way to vectorize search????  regexp search?
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Initialize some arrays and variables
%------------------------------------------------------------------------

% check if strings_to_search is a cell
if ~iscell(strings_to_search)
	% if not, is it a string?
	if ischar(strings_to_search)
		% if so, store it in a cell (size = (1, 1))
		tmpstr = strings_to_search;
		clear strings_to_search;
		strings_to_search = {tmpstr};
		clear tmpstr
	else
		% error!
		error('%s: input strings are in unexpected format', mfilename);
	end
end

%------------------------------------------------------------------------
% find unique values
%------------------------------------------------------------------------
[uniqueText, firstind, locs] = unique(strings_to_search, 'stable');
% uniqueText is a list (cell array) of unique strings
% firstind is a list of indices in strings_search where the unique strings
%	are first found
% locs is a list of where the values in uniqueText occur, 
%	coded by 1, 2, 3, .... nunique

% reorganize as a row vector (for compat. with old findUniqueText)
uniqueText = uniqueText';
% # of unique values found.
NuniqueText = length(uniqueText);
% pre-allocate uniqueIndices
uniqueIndices = cell(NuniqueText, 1);

% loop through # of unique strings
for n = 1:NuniqueText
	% find where each of string index n is found
	uniqueIndices{n} = find(n == locs)';
end
	



%{
% pre-allocate things
Nmarkers = length(strings_to_search);
NuniqueText = 0;
uniqueText = {};
uniqueIndices = {};
textIndex = 1;
testIndex = 1;
runFlag = 1;

% Make a copy of the strings_to_search cell array, stored in testvals to search
% testvals will be used as temporary storage of yet-unmatched strings
testvals = strings_to_search;

%-----------------------------------------------------------------------------
% loop through markers or until runFlag is set to 0 by check within the loop
%-----------------------------------------------------------------------------
while runFlag && (NuniqueText <= Nmarkers)
	% grab new string value to search for
	searchstr = testvals{testIndex};
	% look for matches in original strings_to_search array
	matchindices = strcmp(searchstr, strings_to_search);
	
	% check if any matches were found
	if ~sum(matchindices)
		% no matches found 
		disp(['no match found for ' searchstr])
		% increment testIndex
		testIndex = testIndex + 1;
		% check if testIndex is out of bounds
		if testIndex > length(testvals)
			% stop loop
			runFlag = 0;
		end
		
	else
		% match of searchstr has been found in testvals{}

		% check if new value found by searching for string in uniqueText
		% cell array (which is where the unique strings are stored)
		uniqueCheck = find(strcmp(searchstr, uniqueText)); %#ok<EFIND>
		
		if ~isempty(uniqueCheck)
			% string is not new, so skip ahead to a new test value
			testIndex = testIndex + 1;
			% if new testIndex value is out of bounds, stop loop
			if testIndex > length(testvals)
				runFlag = 0;
			end
		else
			% New (unique) value found
			
			% increment NuniqueText counter
			NuniqueText = NuniqueText + 1;
			% store the current search string in uniqueText cell array
			uniqueText{NuniqueText, 1} = searchstr; %#ok<AGROW>
			% save strings_to_search indices for this search string
			uniqueIndices{NuniqueText, 1} = find(matchindices); %#ok<AGROW>
		
			% rebuild the testvals cell array by removing the instances of the
			% current search string
			newIndices = find(~matchindices);
			if isempty(newIndices)
				% if no new strings were found, stop search
				runFlag = 0;
			else
				% otherwise, assign unique strings to testvals
				testvals = cell(length(newIndices), 1);
				for n = 1:length(newIndices)
					testvals{n} = strings_to_search{newIndices(n)};
				end
				testIndex = 1;
			end
			
		end	% end of "if ~isempty(uniqueCheck)"
	end	% end of "if ~sum(matchindices)"
end	% end of "while runFlag && (NuniqueText <= Nmarkers)"
%}