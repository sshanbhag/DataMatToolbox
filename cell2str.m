function stringOut = cell2str(C, strdelim)
%-----------------------------------------------------------------------------
% stringOut = cell2str(C, strdelim)
%-----------------------------------------------------------------------------
% 						
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neoucom.edu
%------------------------------------------------------------------------
% Created: 21 June, 2011 (SJS)
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Initialize some arrays and variables
%------------------------------------------------------------------------

[Nrows, Ncols] = size(C);
stringOut = cell(Nrows, 1);

if nargin == 2
	switch(upper(strdelim))
		case 'TAB'
			strdelim = '\t';
		case 'SPACE'
			strdelim = ' ';
		case 'COMMA'
			strdelim = ',';
		case 'NONE'
			strdelim = '';
		otherwise
			disp([mfilename ': using default delimiter (space)'])
			strdelim = ' ';
	end
else
	strdelim = ' ';
end


for n = 1:Nrows
	stringOut{n} = '';
	
	for c = 1:Ncols
		
		if isempty(C{n, c})
			tmpstr = '';
			
		elseif ischar(C{n, c})
			tmpstr = sprintf(['%s' strdelim], C{n, c});
			
		elseif isnumeric(C{n,c})
			tmpstr = sprintf(['%f' strdelim], C{n, c});
			
		else
			warnstr = [mfilename ' warning: unknown type for '];
			warnstr = [warnstr sprintf('C{%d. %d} = \n\n', n, c)];
			disp(warnstr);
			disp(C{n, c});
			tmpstr = '?';
		end
		
		stringOut{n} = [stringOut{n} tmpstr];
	end
end

			
			