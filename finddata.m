function varargout = finddata(searchfields, searchstr, ...
															D, fieldnames, varargin)
%------------------------------------------------------------------------
% matchrows = finddata(searchfields, searchstr, D, fieldnames)
% [matchrows, matchdata, allmatches] = finddata(searchfields, searchstr, D, fieldnames)
%------------------------------------------------------------------------
% DataMatToolbox
%------------------------------------------------------------------------
% given searchfields (as single string), finddata will look for the 
% appropriate column in D (using information in fieldnames), and then search
% that column for values that match searchstr.  The parts of D that 
% have matches will be returned in matchdata and the logical vector of 
% matches will be returned in matchrows.
%
% If searchfields and fieldnames are cell arrays, then the individual
% values of each search will be combined using AND - this allows searching
% for multiple congruent values
%
% If a different operation from the default AND is desired, provide the 
% operation as an additional option input:
% 		'CombineMode', 'OR'
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	searchfields		name of column of data in D to search
% 										- OR -
% 							cell list of strings that contain field names 
% 							to search in D
% 							
%	searchstr			string to search for in D
% 										- OR	-
% 							cell list of strings to search in D
% 							
%	D						database cell array (e.g., FullData from FullData6.mat)
%	fieldnames			names of columns in D (usually header.fields from
% 							FullData6.mat file)
% 
%	Optional Input:
% 		'CombineMode'		'OR'		uses OR to combine searches
% 								'AND'		(default) uses AND to combine searches
% 
% Output Arguments:
% 	matchdata			N x (# columns in D) array of data from D that match
% 							searchfields
% 	matchrows			N matches list of rows that match cases in searchstr
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 20 March, 2013 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% check # of strings to search
%------------------------------------------------------------------------
if iscell(searchfields)
	Nsearches = length(searchfields);
else
	% only 1 string/search to use, convert to cell for consistency
	Nsearches = 1;
	searchfields = {searchfields};
	searchstr = {searchstr};
end

%------------------------------------------------------------------------
% Check input varargs
%------------------------------------------------------------------------
if isempty(varargin)
	CombineMode = 'AND';
else
	CombineMode = upper(varargin{1});
end


%------------------------------------------------------------------------
% First, find which column of D contains each searchfield
%------------------------------------------------------------------------
col = finddatacolumn(searchfields, fieldnames);

%------------------------------------------------------------------------
% pre-allocate tmpmatch matrix (logical)
%------------------------------------------------------------------------
[drows, dcols] = size(D);
tmpmatch = false(drows, Nsearches);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
for c = 1:Nsearches
	% look for each string
	% syntax note: by using normal parenthes (instead of
	% {} typically used for cell arrays), matlab returns a cell vector of
	% values.  If {} were used, matlab would return multiple, individual 
	% values and an error would occur using strcmpi
	tmpmatch(:, c) = strcmpi(searchstr{c}, D(:, col(c)));
end

%------------------------------------------------------------------------
% AND or OR the tmpmatch values if Nsearches > 1
%------------------------------------------------------------------------
if Nsearches > 1
	if strcmpi(CombineMode, 'AND')
		matchrows = all(tmpmatch, 2);
	elseif strcmpi(CombineMode, 'OR')
		matchrows = any(tmpmatch, 2);
	else
		error('%s: unknown CombineMode %s', mfilename, CombineMode);
	end
else
	matchrows = tmpmatch;
end

%------------------------------------------------------------------------
% assign outputs.
%------------------------------------------------------------------------
% return rows for matches
if any(nargout == [0 1 2 3 4])
	varargout{1} = matchrows;
end

% return the records that were found
if any(nargout == [2 3 4])
	varargout{2} = D(matchrows, :);
end

% return the raw matches for each searchstr
if any(nargout == [3 4])
	varargout{3} = tmpmatch;
end

% return columns for data
if any(nargout == 4)
	varargout{4} = col;
end

