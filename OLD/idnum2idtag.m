function idtag = idnum2idtag(idnum, cellstr)
%------------------------------------------------------------------------
%idtag = idnum2idtag(idnum, cellstr)
%------------------------------------------------------------------------
% given input cell array of strings cellstr (usually N X 1, 
% but N X M will work), and idnumber code, returns id tag
% corresponding to idnum
%------------------------------------------------------------------------
% Input Arguments:
%	idnum		number
%	cellstr	list of strings to search
%
% Output Arguments:
%	idtag		string
%
%------------------------------------------------------------------------
% See also:  
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 14 June, 2011 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

if (idnum > 0) && (idnum <= length(cellstr))
	idtag = cellstr{idnum};
else
	idtag = -1;
end
