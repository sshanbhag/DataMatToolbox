function [Hvalue, Hmean, HProb] = Hstat(Ddistance)

%Number of random matrixes
RandSets = 10;
%Calculates the significane for a confusion matrix
Nset = Ddistance;
RespNum = length(Nset(:,1));
StimNum = length(Nset(1,:));
Ntot = sum(sum(Nset));
Nlength = RespNum * StimNum;

TotalSum = 0;
% Measure the information content of the stimulus matrix;
for i = 1:RespNum,
	for ii = 1:StimNum,
		StimSum = 0;
		RespSum = 0;
		TempSum = 0;
		for a = 1:StimNum
			StimSum = StimSum + Nset(a,ii);
		end;
		for b = 1:RespNum
			RespSum = RespSum + Nset(i,b);
		end;
		if Nset(i,ii) == 0;
			TempSum =0;
		elseif and(StimSum ==0,RespSum~=0);
			TempSum = Nset(i,ii)*(log(Nset(i,ii))-log(RespSum) + log(Ntot));
		elseif and(StimSum ~=0,RespSum==0);
			TempSum = Nset(i,ii)*(log(Nset(i,ii))-log(StimSum)+ log(Ntot));
		else
			TempSum = Nset(i,ii)*(log(Nset(i,ii))-log(StimSum)-log(RespSum) + log(Ntot));
		end;
		TotalSum = TempSum + TotalSum;
	end;
end;
Hvalue = (1/Ntot)*TotalSum;

if nargout == 1
	return
end

%Creates some random matrixes to compare to the real one;
Index = zeros(RandSets, Nlength);
for i = 1:RandSets;
    Index(i,:) = randperm(Nlength);
end;
for r = 1:RandSets,
    Ntemp=  reshape(Nset,1,Nlength);
    Ntemp = Ntemp(1,Index(r,:));
    Nset2 = reshape(Ntemp,StimNum,RespNum);
    TotalSum2 = 0;
	for i = 1:StimNum,
		 for ii = 1:RespNum,
			  StimSum = 0;
			  RespSum = 0;
			  TempSum = 0;
			  for a = 1:StimNum
					StimSum = StimSum + Nset2(a,ii);
			  end;
			  for b = 1:RespNum
					RespSum = RespSum + Nset2(i,b);
			  end;
			  if Nset2(i,ii) == 0;
					TempSum =0;
			  elseif and(StimSum ==0,RespSum~=0);
					TempSum = Nset2(i,ii)*(log(Nset2(i,ii))-log(RespSum) + log(Ntot));
			  elseif and(StimSum ~=0,RespSum==0);
					TempSum = Nset2(i,ii)*(log(Nset2(i,ii))-log(StimSum)+ log(Ntot));
			  else
					TempSum = Nset2(i,ii)*(log(Nset2(i,ii))-log(StimSum)-log(RespSum) + log(Ntot));
			  end;

			  TotalSum2 = TempSum + TotalSum2;
		 end;
		 HRand(:,r) = (1/Ntot)*TotalSum2;
	end;
end;
%This give the confidence intervals of the H value from the bootsrap data

HSd = std(HRand);
Hmean = mean(HRand);
Sigvalue = [{'<0.01'},{'<0.10'}, {'>0.30'}, {'<0.05'},{'0.01'}];
HSe(1) = Hmean-HSd*3;
HSe(2) = Hmean-HSd*2;
HSe(3) = Hmean;
HSe(4) = Hmean + HSd*2;
HSe(5) = Hmean + HSd*3;
if Hvalue <= HSe(2),
    HProb = [0.05];
elseif Hvalue >= HSe(4),
    HProb = [0.05];
else
    HProb = [0.2];
end;
disp(HProb);