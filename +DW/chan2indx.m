function indx = chan2indx(chanstr)
if strcmpi(chanstr, 'L')
	indx = 1;
else
	indx = 2;
end
