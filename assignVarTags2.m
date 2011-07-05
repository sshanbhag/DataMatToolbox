
function Var = assignVarTags(Stimstruct, VarTags)
%------------------------------------------------------------------------
% Var = assignVarTags(Stimstruct, VarTags)
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%
%
%------------------------------------------------------------------------
% See:  
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 4 July, 2011 (SJS)
% 	- uses code snipped from buildStimulusStruct.m
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% get the number of tags to assign
ntags = length(VarTags);

Var = struct('name', [], 'values', [], 'indices', []);


%%%%%%%
%{
need different way in order to handle character and numeric values.
at present, strings get mangled

Also, need to change algorithm - create list of tags and then search for unique values
%}





% loop through the tags
varindex = 0;
for n = 1:ntags

	clear tmpvar

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
				tmpind = find(tmp == uniqtmp(u));
				tmpvar(1).indices = Stimstruct.Indices(tmpind);
			end
		end

		% get value for Left channel, store in tmpvar(2)
		tmp = Stimstruct.([VarTags{n} 'L']);
		uniqtmp = unique(tmp);
		if ~isempty(uniqtmp)
			tmpvar(2).name = [VarTags{n} 'L'];
			tmpvar(2).values = uniqtmp;

			for u = 1:length(uniqtmp)
				tmpind = find(tmp == uniqtmp(u));
				tmpvar(2).indices = Stimstruct.Indices(tmpind);
			end
		end

	else
		% only 1 channel used
		tmp = Stimstruct.([VarTags{n} Stimstruct.Channel]);
		uniqtmp = unique(tmp);
		if ~isempty(uniqtmp)
			tmpvar.name = [VarTags{n} Stimstruct.Channel];
			tmpvar.values = uniqtmp;
			for u = 1:length(uniqtmp)
				tmpind = find(tmp == uniqtmp(u));
				tmpvar.indices = Stimstruct.Indices(tmpind);
			end
		end
	end

	% assign tmpvar to output struct
	if exist('tmpvar')
		varindex = varindex + 1;
		Var(varindex) = tmpvar;
	end

end






