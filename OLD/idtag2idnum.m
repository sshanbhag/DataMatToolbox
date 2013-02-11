function [idnum, idvect] = idtag2idnum(idtag, cellstr)
%------------------------------------------------------------------------
%[idnum, idvect] = idtag2idnum(idtag, cellstr)
%------------------------------------------------------------------------
% given input cell array of strings cellstr (usually N X 1, 
% but N X M will work), returns numeric id location idnum of search string 
% idtag
%------------------------------------------------------------------------
% Input Arguments:
%	idtag			string to match in cellstr
%	cellstr		cell vector (or array) of strings
%
% Output Arguments:
%	idnum			linear location of idtag in cellstr
%	idvect		vector of size(cellstr) indicating matches for idtag in cellstr
%
%------------------------------------------------------------------------
% See also:  idnum2idtag
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

idvect = strcmp(cellstr, idtag);
idnum = find(idvect);