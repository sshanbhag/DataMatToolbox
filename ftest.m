function varargout = ftest(varargin)

fprintf('isempty(varargin) = %d\n', isempty(varargin));
Nargs = length(varargin);
fprintf('length(varargin) = %d\n', Nargs);
for n = 1:Nargs
	fprintf('varargin{%d} = %s\n', n, inputname(n))
end
if nargout
	varargout{1} = Nargs;
end
