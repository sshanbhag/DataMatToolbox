%% load matfile (if FullData isn't already in workspace)
if ~exist('FullData', 'var')
	load('/Users/sshanbhag/Work/Data/DataWave/batmat/FullData6.mat');
end

%% some ways to parse data


%% by File Type

% First, find which column is File Type
%	One way to do this is by searching for the string 'File Type'
%	in the header field names struct - strcmpi can do this and will return
%	a vector that has 1 where there is a match and 0 for non-matching strings
%	The location (index) of matches is then determined using the find() function
col = find(strcmpi('File Type', header.fields))

% if not found, throw an error
if isempty(col)
	error('string "File Type" not found in header.fields!')
end

% next, define a string to look for
searchstr = 'strings block_Sorted_255.mat';

% look for this string - use the header.isnumber value to determine how to
% search

if header.isnumber(col) == 1
	% values are a number`
	
else
	% values are strings
	% look for matches.  
	% syntax note: by using normal parenthes (instead of
	% {} typically used for cell arrays), matlab returns a cell vector of
	% values.  If {} were used, matlab would return multiple, individual 
	% values and an error would occur using strcmpi
	matchrows = strcmpi(searchstr, FullData(:, col));
	% check if any were found
	if isempty(matchrows)
		error('No matches for string %s found in FullData column %d!', ...
																				searchstr, col);
	end
	
	% store the records that were found
	MatchData = FullData(matchrows, :);
end


%% Look for 2 conditions

% We want only complete, responsive units.
% find column for the DeleteRow data
deletecol = find(strcmpi('DeleteRow', header.fields))
% find column for the Unresponsive Unit data
unrespcol = find(strcmpi('Unresponsive Unit', header.fields))

% if not found, throw an error
if any( [ isempty(deletecol) isempty(unrespcol) ] )
	error('column string not found in header.fields!')
end


% look for matches to keep from DeleteRow ('0' means a complete test)  
deleterows = strcmpi('0', FullData(:, deletecol));
% look for matches to keep from Unresponsive Unit (keep responsive units)
unresprows = strcmpi('0', FullData(:, unrespcol));

%  use logical AND operator (&) to look for rows for which both conditions
% are true
rows_to_keep = deleterows & unresprows;

MatchData = FullData(rows_to_keep, :);

%% and combine the complete responsive units with Type 1 data

% look for column
typecol = find(strcmpi('Type', header.fields));
% look for rows
typerows = strcmpi('1', FullData(:, typecol));
% AND the three conditions
rows3 = deleterows & unresprows & typerows;

% get data
Complete1 = FullData(rows3, :);

%% get unique unit numbers
AllUnitNumbers = zeros(length(Complete1), 1);
for n = 1:length(Complete1)
	AllUnitNumbers(n) = str2num(Complete1{n, 2});
end

% find unique unit numbers
UnitNumbers = unique(AllUnitNumbers);

%% get responsiveness 
attenvals = {'0', '20', '40'};
for n = 1:length(UnitNumbers)
	% find Complete1 data for current unit
	unitrows = find(strcmpi(num2str(UnitNumbers(n)), Complete1(:, 2)))
	
	UnitAttenResp{n} = cell(3, 1);
	
	
	
	
	return
end



