function Var = assignVarTags(Stimstruct, VarTags)
%------------------------------------------------------------------------
% Var = assignVarTags(Stimstruct, VarTags)
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Stimstruct			DataWave toolbox stimulus structure
%	VarTags				cell string list of tags to search in Stimstruct
% 
% Output Arguments:
%	Var					Variable variation structure array
% 								length is equivalent to length of VarTags elements
% 		name
% 		values
% 		indices
%
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 4 July, 2011 (SJS)
% 	- uses code snipped from buildStimulusStruct.m
%
% Revisions:
%	6 July, 2011 (SJS): 
% 		-	updated documentation
% 		-	made changes to allow for different types of search vars
% 			e.g., strings, numeric vectors, cell vectors
%	14 July, 2011 (SJS): some tweaks for speed
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% get the number of tags to assign
ntags = length(VarTags);

Var = struct('name', [], 'values', [], 'indices', []);

% loop through the tags
varindex = 0;
for n = 1:ntags

	clear tmpvar;
	clear uniqmp uniqindx nuniq;


	% check if both channels are used and assign variable tags and values
	% accordingly
	if Stimstruct.Channel == 'B'
		% both channels used

		% get the value for the current tag, right channel
		tmp = Stimstruct.([VarTags{n} 'R']);
		% look for unique values
		uniqtmp = unique(tmp);
		if ~isempty(uniqtmp)
			% if unique, store it in 
			tmpvar(1).name = [VarTags{n} 'R'];
			tmpvar(1).values = uniqtmp;

			for u = 1:length(uniqtmp)
% 				tmpind = find(tmp == uniqtmp(u));
				tmpvar(1).indices{u} = Stimstruct.Indices(tmp == uniqtmp(u));
			end
		end

		% get value for Left channel, store in tmpvar(2)
		tmp = Stimstruct.([VarTags{n} 'L']);
		uniqtmp = unique(tmp);
		if ~isempty(uniqtmp)
			tmpvar(2).name = [VarTags{n} 'L'];
			tmpvar(2).values = uniqtmp;

			for u = 1:length(uniqtmp)
% 				tmpind = find(tmp == uniqtmp(u));
				tmpvar(2).indices{u} = Stimstruct.Indices(tmp == uniqtmp(u));
			end
		end

	else
		% only 1 channel used, only get info for that channel
		tmp = Stimstruct.([VarTags{n} Stimstruct.Channel]);
		
		% need to handle tags differently depending on data types
		if ischar(tmp)
			[uniqtmp, uniqindx, nuniq] = findUniqueText(tmp);
		elseif iscell(tmp)
			if ischar(tmp{1})
				[uniqtmp, uniqindx, nuniq] = findUniqueText(tmp);
			else
				[uniqtmp, uniqindx, nuniq] = findUniqueCellRows(tmp);
			end
		else
			[uniqtmp, uniqindx, nuniq] = findUniqueValues(tmp);
		end
		
		if exist('uniqtmp', 'var')
			if ~isempty(uniqtmp)
				tmpvar.name = [VarTags{n} Stimstruct.Channel];
				tmpvar.values = uniqtmp;
				tmpvar.indices = uniqindx;
			end
		end
	end

	% assign tmpvar to output struct
	if exist('tmpvar', 'var')
		varindex = varindex + 1;
		Var(varindex) = tmpvar;
	end

end






