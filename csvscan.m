function [out, nread] = csvscan(in, multdelim)
	tmp = textscan(	in, '%s', ...
							'Delimiter', ',', ...
							'MultipleDelimsAsOne', multdelim);
	out = tmp{1};
	nread = length(out);
	clear tmp
