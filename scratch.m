

%{
% plot, for each unit, the information -vs- cost
for u = 1:length(Ddata)
	h = Ddata(u).Hvalue;
	figure(u)
	plot(CostVals, h, '.-')
	xlabel('cost/sec')
	ylabel('H (bits)')
	[a, aind] = max(h);
	hold on
		plot(CostVals(aind), a, 'ro')
	hold off
	title({	sprintf('File %s', filename1),	...
				sprintf('Unit %d', u + 1),			...
				sprintf('qmax=%.2f', CostVals(aind))	}, ...
				'Interpreter', 'none');
end

% store Ddata in OutData cell array
OutData{findex} = Ddata;
clear Ddata;
	


%}

