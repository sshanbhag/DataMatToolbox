function cols = finddatacolumn(searchfields, fieldnames)
%------------------------------------------------------------------------
% cols = finddatacolumn(searchfields, fieldnames)
%------------------------------------------------------------------------
% DataMatToolbox
%------------------------------------------------------------------------
% given searchfield (as single string), finddatacolumn will look for the 
% appropriate column in D (using information in fieldnames), and return the
% column number
%------------------------------------------------------------------------
% Input Arguments:
% 	searchfields		name of column of data in D to search
% 										- OR -
% 							cell list of strings that contain field names 
% 							to search in D
% 							
%	fieldnames			names of columns in D (usually header.fields from
% 							FullData6.mat file)
% 										- OR -
% 							can simply pass in the header data struct
% 
% Output Arguments:
%	cols						list of column numbers
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
%	12 Mar 2015 (SJS): ability to use struct as fieldnames input added
%							to help info
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
end

%------------------------------------------------------------------------
% check if fieldnames is struct - if so, extract fieldnames
%------------------------------------------------------------------------
if isstruct(fieldnames)
	tmpf = fieldnames.fields;
	clear fieldnames;
	fieldnames = tmpf;
	clear tmpf;
end

%------------------------------------------------------------------------
% First, find which column of D contains each searchfield
%------------------------------------------------------------------------
%	One way to do this is by searching for the string 'File Type'
%	in the header field names struct - strcmpi can do this and will return
%	a vector that has 1 where there is a match and 0 for non-matching strings
%	The location (index) of matches is then determined using the find() function
%------------------------------------------------------------------------

cols = zeros(Nsearches, 1);
for c = 1:Nsearches
	tmpcol = find(strcmpi(searchfields{c}, fieldnames));
	% if not found, throw an error
	if isempty(tmpcol)
		error('%s: string "%s" not found in fieldnames!', ...
																	mfilename, searchfields{c});
	else
		cols(c) = tmpcol;
	end
end




