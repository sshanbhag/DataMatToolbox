
nFiles = length(BBN);
for n = 1:nFiles
	nUnitList(n) = BBN{n}.D.Info.Nunits;
	blist{n} = BBN{n}.D.Info.UnitList;
end

maxNunits = max(nUnitList);
minNunits = min(nUnitList);

UnitList = zeros(min(nUnitList, 2);
UnitIndex = UnitList;

u = 0;
for p = 1:maxNunits
	lFlag = 0;
	for l = 1:LFH.D.Info.Nunits
		% check if probe # (col 1) and cluster # (col 2) match
		if all(blist(p, 1:2) == llist(l, 1:2))
			lFlag = l;
		end
	end
	% if match found, store the unit information
	if lFlag
		u = u + 1;
		UnitList(u, :) = blist(p, 1:2);
		UnitIndex(u, :) = [p lFlag];
	end
end


return

%------------------------------------------------------------------------
% make sure # of units match
%------------------------------------------------------------------------
%	store matching units information in UnitList
%------------------------------------------------------------------------
if BBN.D.Info.Nunits ~= LFH.D.Info.Nunits
	warning('%s: # Units mismatch!!!!!!!!!!!\n', mfilename)
	fprintf('\tBBN file has %d\tunits\n', BBN.D.Info.Nunits);
	fprintf('\tLFH file has %d\tunits\n\n', LFH.D.Info.Nunits);

	% make some local copies of the UnitList struct from each file
	blist = BBN.D.Info.UnitList;
	fprintf('\tBBN Unit Info:\n')
	fprintf('\t\tUnit\t\tProbe\t\tCluster\n')
	for p = 1:BBN.D.Info.Nunits
		fprintf('\t\t%d\t\t\t%d\t\t\t%d\n', blist(p, 3), blist(p, 1), blist(p, 2) );
	end

	llist = LFH.D.Info.UnitList;
	fprintf('\tLFH Unit Info:\n')
	fprintf('\t\tUnit\t\tProbe\t\tCluster\n')
	for p = 1:LFH.D.Info.Nunits
		fprintf('\t\t%d\t\t\t%d\t\t\t%d\n', llist(p, 3), llist(p, 1), llist(p, 2) );
	end

	% if mismatch # of units, find units that are common to both
	% i.e., same probe and cluster #
	% store this list in UnitList which is a n X 3 matrix with
	%  [ Probe#		Cluster#	] as columns
	% and store unit numbers in 
	% UnitIndex[BBN.D.Info.UnitList(:, 3) LFH.D.Info.UnitList(:, 3)]
	if BBN.D.Info.Nunits > LFH.D.Info.Nunits
		% total # of matching units MUST be the lower of the 
		% # of units in the two files... in this case, LFH
		UnitList = zeros(LFH.D.Info.Nunits, 2);
		UnitIndex = UnitList;
		u = 0;
		for p = 1:BBN.D.Info.Nunits
			lFlag = 0;
			for l = 1:LFH.D.Info.Nunits
				% check if probe # (col 1) and cluster # (col 2) match
				if all(blist(p, 1:2) == llist(l, 1:2))
					lFlag = l;
				end
			end
			% if match found, store the unit information
			if lFlag
				u = u + 1;
				UnitList(u, :) = blist(p, 1:2);
				UnitIndex(u, :) = [p lFlag];
			end
		end
	else
		% fewer units in BBN file, so that is the # of matching units
		UnitList = zeros(BBN.D.Info.Nunits, 2);
		UnitIndex = UnitList;
		u = 0;
		for p = 1:LFH.D.Info.Nunits
			bFlag = 0;
			for l = 1:BBN.D.Info.Nunits
				if all(llist(p, 1:2) == blist(l, 1:2))
					bFlag = l;
				end
			end

			if bFlag
				u = u + 1;
				UnitList(u, :) = llist(p, 1:2);
				UnitIndex(u, :) = [bFlag p];
			end
		end
	end

	Nunits = length(UnitList(:, 1));
else
	% units match so set UnitList to one of the lists from the 
	% D structures (for no particular reason, BBN is used)
	Nunits = BBN.D.Info.Nunits;
	UnitList = BBN.D.Info.UnitList(:, 1:2);
	UnitIndex = [BBN.D.Info.UnitList(:, 3) LFH.D.Info.UnitList(:, 3)];
end
