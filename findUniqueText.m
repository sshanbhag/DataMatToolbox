function [uniqueText, uniqueIndices, NuniqueText] = findUniqueText(strings_to_search)
%-----------------------------------------------------------------------------
%[uniqueText, uniqueIndices, NuniqueText] = findUniqueText(strings_to_search)
%-----------------------------------------------------------------------------
% DataMatToolbox
%-----------------------------------------------------------------------------
% Finds unique strings within the list of strings_to_search
% 
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
% 	 -	revised algorithm to use unique() MATLAB function.  this is
% 		faster for longer lists!
%	21 Mar 2013 (SJS): cleaned up code
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
% need to account for differences in unique() function for matlab versions
mver = version('-release');
if strcmp(mver(1:4), '2012')
	[uniqueText, firstind, locs] = unique(strings_to_search, 'first');
% 	unique2012.strings = strings_to_search;
% 	unique2012.uniqueText = uniqueText;
% 	unique2012.firstind = firstind;
% 	unique2012.locs = locs;
% 	save unique2012.mat unique2012
else
	[uniqueText, firstind, locs] = unique(strings_to_search, 'first');
% 	unique2010b.strings = strings_to_search;
% 	unique2010b.uniqueText = uniqueText;
% 	unique2010b.firstind = firstind;
% 	unique2010b.locs = locs;
% 	save unique2010b.mat unique2010b
end

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

