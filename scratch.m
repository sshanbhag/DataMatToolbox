for d = 1:Ndata
	validDataIndices = validBBNList{d, 1};
	validUnitIndices = validBBNList{d, 5};
	validAttenIndices = validBBNList{d, 7};
	[NvalidUnits, tmp] = size(validUnitIndices);
	
	
	
	for u = 1:NvalidUnits
		bbnData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitList
		for c = 1:Nconditions
			clustn(u, c) = bbnData{validDataIndices(c)}.UnitData(validUnitIndices(u, c)).UnitInfo.cluster;
		end
	end
	
	clustn
	
return
	
end