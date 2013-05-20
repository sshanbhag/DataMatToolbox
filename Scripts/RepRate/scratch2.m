% clear everything
clear all

% create some variables
a_voltage = [1 2 3];
b_voltage = [2 34 4 5];
b_current = [4 5 64 56 456];

% get list of variables
varlist = who


% loop through varlist
for n = 1:length(varlist)
	testval = strfind(varlist{n}, '_voltage')
	
	% check if something matched
	if ~isempty(testval)
		% if so, display square of the values
		eval(sprintf('%s.^2', varlist{n}))
		eval(sprintf('figure; plot(%s)', varlist{n}))
		
	end
	
end