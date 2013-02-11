function fileList = find_condition_instances(stringList, fullList)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% find matching filenames/conditions
%
%	store results in bbnList
%		Col 1: match_index (index into bbn_idstr, and files)
%		Col 2: 
%--------------------------------------------------------------------------
tmpList = stringList;
match_index = 0;
runFlag = 1;
n = 1;
while runFlag
	% find all instances of tmpList{n} in the list tmpList
	tmpsrch = strcmp(tmpList{n}, tmpList);

	% if more than 1 are found, store info
	if sum(tmpsrch) > 1
		% a match was found (apart from with itself...), so increment
		% match_index counter
		match_index = match_index + 1;
		% redo instance search, but this time in the full idstr list
		% (not the truncated one in tmpsrch)
		srchout = strcmp(tmpList{n}, stringList);
		% store data in fileList
		
		% store indices into fullList
		fileList{match_index, 1} = find(srchout);
		% store matched strings
		fileList{match_index, 2} = stringList(srchout);
		% store matched full list members (corresponding to matched strings)
		fileList{match_index, 3} = fullList(srchout);
		
		% refresh tmpidstr to eliminate instances of tmpidstr{n}
		% in the running list - this will eliminate duplicate
		% searches.
		tmpList = tmpList(~tmpsrch);
		% reset n to 1
		n = 1;
	else
		% no match was found
		n = n + 1;
	end

	% trap errors/end

	if isempty(tmpList)
		runFlag = 0;
	elseif length(tmpList) == 1
		runFlag = 0;
	elseif n > length(tmpList)
		runFlag = 0;
	else
		runFlag = 1;
	end

end