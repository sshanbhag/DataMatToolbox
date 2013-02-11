function [M, L] = match_atten(Cdata)

Ncond = length(Cdata);

if Ncond == 1
	warning('DATAMAT:redundant', '%s: no need to search, only 1 condition found!', mfilename)
	return
end

nAttenList = zeros(Ncond, 1);
for c = 1:Ncond
	nAttenList(c) = length(Cdata{c});
end

% find max, min of # units
[maxNatten, maxIndex] = max(nAttenList);

% loop through the "base" array that has the most atten vals

% get the "base" array, and build list of test arrays
baseC = Cdata{maxIndex};
testindex = 1:Ncond;
testindex = testindex(testindex ~= maxIndex);

ntest = length(testindex);
testC = cell(ntest, 1);
for n = 1:ntest
	testC{n} = Cdata{testindex(n)};
end

M = zeros(maxNatten, Ncond);

% loop through max # units (base array)
for u = 1:maxNatten
	M(u, maxIndex) = u;
	% loop through test arrays
	for t = 1:ntest
		testval = false(nAttenList(testindex(t)), 1);
		testout = zeros(nAttenList(testindex(t)), 1);
		
		% loop through units for current test array
		for r = 1:nAttenList(testindex(t))	
			% compare base pair to test pair
			testval(r) = isequal(baseC(u), testC{t}(r));
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
			warning('DATAMAT:unitnum', '%s: multiple identical attenuation values found!', mfilename)
			M(u, testindex(t)) = testout(matchind(1));
		elseif ~isempty(matchind)
			M(u, testindex(t)) = testout(matchind);
		end	
	end
end

allM = all(M, 2);

if ~any(allM)
	warning('%s: no common attenuation values', mfilename);
	fprintf('\tUsing lowest value from each file\n')
	M = ones(1, Ncond);
	L = ones(1, Ncond);
else
	[nr, nc] = size(M);
	u = 0;
	for n = 1:nr
		if allM(n)
			u = u + 1;
			L(u, :) = M(n, :);
		end
	end
end
	



