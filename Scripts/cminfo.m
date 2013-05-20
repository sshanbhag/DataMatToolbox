function [H, Hs, Hr, Hsr] = cminfo(CM)

[nrows, Nstimuli] = size(CM);
if nrows ~= Nstimuli
	error('%s: CM must equal # of rows and columns', mfilename);
end

respcounts = round(sum(CM, 2));
Nresp = respcounts(1);
if ~all(Nresp == respcounts)
	display(respcounts)
	error('%s: sum of all rows must be equal', mfilename);
end

% marginal prob for stimuli
Ps = ones(Nstimuli, 1) ./ Nstimuli;

% marginal for responses
% sum along columns
rcount = sum(CM, 1);
% divide by grand total
Pr = rcount ./ sum(sum(CM));

% joint prob
Psr = (CM ./ Nresp) ./ Nstimuli;

% now compute entropies

Psl = log2(Ps);
Psl(Psl == -Inf) = 0;
Psl(Psl == Inf) = 0;
Hs = -1 * sum( Ps .* Psl);

Prl = log2(Pr);
Prl(Prl == -Inf) = 0;
Prl(Prl == Inf) = 0;
Hr = -1 * sum( Pr .* Prl);

Psrl = log2(Psr);
Psrl(Psrl == -Inf) = 0;
Psrl(Psrl == Inf) = 0;
Hsr = -1 * sum(sum( Psr .* Psrl));

% total entropy (MI)
H = Hs + Hr - Hsr;


