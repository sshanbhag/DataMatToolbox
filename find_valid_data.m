function validList = find_valid_data(Data, List)
%--------------------------------------------------------------------------
% validList = find_valid_data(Data, List)
%--------------------------------------------------------------------------
% limit list to analyze to those files for which there are all three
% conditions (1, 2, 3)
%--------------------------------------------------------------------------
% valid___List{}
% 	{n, 1}	->	index into Data and bbnFiles
% 	{n, 1}	->	file search strings (no _BBN or _LFH or condition value)
% 	{n, 2}	->	full .mat file name
% 	{n, 4}	->	condition #
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

n = 0;
for bIndx = 1:length(List)
	clist = List{bIndx, 4};
	if (~isempty(clist)) && (length(clist) >= 3)
		% see if desired conditions are in the list
		ctest = zeros(3, length(clist));
		ctest(1, :) = (clist == 1);
		ctest(2, :) = (clist == 2);
		ctest(3, :) = (clist == 3);
		if sum(sum(ctest)) == 3
			for c = 1:3
				cind(c) = find(ctest(c, :));
			end

			% if so, sort and store the file information
			n = n + 1;
			validList(n, :) = List(bIndx, :);
			for c = 1:3
				for j = 1:4
					validList{n, j} = validList{n, j}(cind);
				end
			end
		end
	end
end


% loop through file list
for fIndx = 1:length(validList)
	% get # of files from column 3 of valid list
	bfiles = validList{fIndx, 3};
	% get # of indices 
	bindex = validList{fIndx, 1};

	if isempty(bfiles)
		error('%s: empty bfiles', mfilename)
	end
	
	Ncond = length(bindex);

	% load TMP file data into list
	TMP = cell(1, Ncond);
	TMPlist = cell(1, Ncond);
	for n = 1:Ncond
		TMP{n} = Data{bindex(n)};
		% columns of blist arrays are probe #, cluster # and unit #
		TMPlist{n} = TMP{n}.Info.UnitList;
	end
	
	[unitM, matches] = match_units(TMPlist);
	if isempty(unitM)
		warning('%s: no unit matches found', mfilename)
	end
	validList{fIndx, 5} = matches;
	validList{fIndx, 6} = unitM;
end
	

%--------------------------------------------------------------------------
% determine number of Attenuation levels by looking at the Stimulus.Var
% information in each condition's data structure
%--------------------------------------------------------------------------
for fIndx = 1:length(validList)
	bfiles = validList{fIndx, 3};
	bindex = validList{fIndx, 1};

	if isempty(bfiles)
		error('%s: empty bfiles', mfilename)
	end
	
	Ncond = length(bindex);

	% load TMP file data into list
	TMP{c} = cell(1, Ncond);
	TMPlist = cell(1, Ncond);
	for c = 1:Ncond
		TMP{c} = Data{bindex(c)};
		% columns of TMPlist arrays are probe #, cluster # and unit #
		TMPlist{c} = TMP{c}.Stimulus.Var;
		% # of attenuation values in TMP file
		nVar = length(TMPlist{c});

		attenindex = 0;
		for n = 1:nVar
			% find R channel attenuation
			if strcmp(TMPlist{c}(n).name, 'AttenuationR')
				attenindex = n;
			end
		end
		
		% if found, store the values
		if attenindex
			attenVals{c} = TMPlist{c}(attenindex).values;
		else
			% if not, error
			error('%s: attenuation not found in Var list for file %s', mfilename, bfiles{c});
		end
		Natten(c) = length(attenVals{c});				
	end

	[m, l] = match_atten(attenVals);
	validList{fIndx, 7} = l;
	validList{fIndx, 8} = m;
end

