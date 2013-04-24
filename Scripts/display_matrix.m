function display_matrix(aMatrix, varargin)

if isempty(varargin)
	fmt = '%f';
else
	fmt = varargin{1};
end


[nr, nc] = size(aMatrix);

fprintf('\t%s = \n', inputname(1));
for r = 1:nr
	fprintf('\t\t');
	pstr = sprintf('fprintf(''%s\t'', aMatrix(%d, :))', fmt, r);
	eval(pstr)
	fprintf('\n');
	
end

	