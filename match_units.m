function [M, L] = match_units(Cdata)

Ncond = length(Cdata);

if Ncond == 1
	warning('DATAMAT:redundant', '%s: no need to search, only 1 condition found!', mfilename)
	return
end

nUnitList = zeros(Ncond, 1);
for c = 1:Ncond
	nUnitList(c) = length(Cdata{c});
end

% find max, min of # units
[maxNunits, maxIndex] = max(nUnitList);

% loop through the "base" array that has the most units
% This will ensure that the max # of possible probe/unit combinations will
% be tested.

% get the "base" array, and build list of test arrays
baseC = Cdata{maxIndex};
testindex = 1:Ncond;
testindex = testindex(testindex ~= maxIndex);

ntest = length(testindex);
testC = cell(ntest, 1);
for n = 1:ntest
	testC{n} = Cdata{testindex(n)};
end

M = zeros(maxNunits, Ncond);

% loop through max # units (base array)
for u = 1:maxNunits
	M(u, maxIndex) = u;
	% loop through test arrays
	for t = 1:ntest
		testval = false(nUnitList(testindex(t)), 1);
		testout = zeros(nUnitList(testindex(t)), 1);
		
		% loop through units for current test array
		for r = 1:nUnitList(testindex(t))	
			% compare base pair to test pair
			testval(r) = isequal(baseC(u, 1:2), testC{t}(r, 1:2));
			%{
			fprintf('u(%d) vs. t(%d)r(%d) -> ', u, t, r);
			fprintf('%d %d\t\t', baseC(u,1), baseC(u,2));
			fprintf('%d %d\t', testC{t}(r,1), testC{t}(r,2));
			fprintf('=\t%d r(%d)\n', testval(r), r);
			%}
			% if found, save the index
			if testval(r)
				testout(r) = r;
			else
				testout(r) = 0;
			end
		end
		
		matchind = find(testout > 0);
		if length(matchind) > 1
			warning('DATAMAT:unitnum', '%s: multiple units found!', mfilename)
			M(u, testindex(t)) = testout(matchind(1));
		elseif ~isempty(matchind)
			M(u, testindex(t)) = testout(matchind);
		end	
	end
end

if nargout > 1
	allM = all(M, 2);
	
	[nr, nc] = size(M);
	u = 0;
	for n = 1:nr
		if allM(n)
			u = u + 1;
			L(u, :) = M(n, :);
		end
	end
end

	



