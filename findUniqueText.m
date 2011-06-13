function [uniqueText, uniqueIndices, NuniqueText] = findUniqueText(strings_to_search)
%-----------------------------------------------------------------------------
%[uniqueText, uniqueIndices, NuniqueText] = findUniqueText(strings_to_search)
%-----------------------------------------------------------------------------
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
% 					- store indices of searchstr in strings_to_search in the uniqueIndices list
% 					- rebuild the testvals list by removing instances of the current
% 					  searchstr
% 						
%-----------------------------------------------------------------------------
% NuniqueText		# of unique strings in strings_to_search
% uniqueText		unique values of strings in strings_to_search, cell array
% uniqueIndices	indices where each of the unique strings are found in strings_to_search
% testvals			temporary storage of unmatched strings
%-----------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 13 June, 2011 (SJS)
% 	- uses code snipped from loadDWStimData.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% Initialize some arrays and variables
Nmarkers = length(strings_to_search);
NuniqueText = 0;
uniqueText = {};
uniqueIndices = {};
textIndex = 1;
testIndex = 1;
runFlag = 1;


% Make a copy of the strings_to_search cell array, stored in testvals to search
testvals = strings_to_search;

% loop through markers or until runFlag is set to 0 by check within the loop
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
		uniqueCheck = find(strcmp(searchstr, uniqueText));
		
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
			uniqueText{NuniqueText, 1} = searchstr;
			% save strings_to_search indices for this search string
			uniqueIndices{NuniqueText, 1} = find(matchindices);
		
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




