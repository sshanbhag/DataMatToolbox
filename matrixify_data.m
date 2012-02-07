function Sout = matrixify_data(Sin)

Nunits = length(Sin);
Nconditions = length(Sin(1).condition);
Nwin = length(Sin(1).Count_mean{1});
Nbgwin = length(Sin(1).BG{1});


tmp = zeros(Nunits, Nconditions);
Sout.unitnum = tmp;
Sout.probe = tmp;
Sout.cluster = tmp;
Sout.BG_total = tmp;
Sout.Net_mean = tmp;
Sout.Net_std = tmp;

tmp = zeros(Nunits, Nwin, Nconditions);
Sout.Count_mean = tmp;
Sout.Count_std = tmp;

Sout.BG_mean = zeros(Nunits, Nbgwin, Nconditions);
Sout.BG_std = zeros(Nunits, Nbgwin, Nconditions);


for b = 1:Nunits
	Sout.unitnum(b, :) = Sin(b).unit;
	Sout.probe(b, :) = Sin(b).probe;
	Sout.cluster(b, :) = Sin(b).cluster;
	Sout.atten(b, :) = Sin(b).AttenVals;
	Sout.BG_total(b, :) = Sin(b).BG_count;
	Sout.Net_mean(b, :) = Sin(b).Net_mean;
	Sout.Net_std(b, :) = Sin(b).Net_std;
	
	for c = 1:Nconditions
		Sout.Count_mean(b, :, c) = Sin(b).Count_mean{c};
		Sout.Count_std(b, :, c) = Sin(b).Count_std{c};
		Sout.BG_mean(b, :, c) = Sin(b).BG_mean{c};
		Sout.BG_std(b, :, c) = Sin(b).BG_std{c};
	end
end
