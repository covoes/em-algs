function [mdl t2] = getMDL(resStruct, r, t, nObj, data)
	t2 = t;
	if isempty(resStruct) || isempty(resStruct{r})
		mdl = NaN;
		return
	end

	if isempty(resStruct{r}{t})
		if t > 1
			%find the last one 
			for t2=t-1:-1:1
				if	~isempty(resStruct{r}{t2})
					mdl = getMDLIntern(resStruct{r}{t2}, nObj, data);
					return;
				end
			end
			mdl = NaN;
		else
			mdl = NaN;
		end
	else
		mdl = getMDLIntern(resStruct{r}{t},nObj, data);
	end
end

function mdl = getMDLIntern(resStruct, nObj, data)
	if isfield(resStruct,'medoids')
		gmmObj =  resultsGAKREMToGMM(resStruct, data);
		%to simplify things (avoiding a estimation step) we use a fake object 
		%and obtain the log likelihood using the cluster function
		gmmObjFake = struct('NComponents', gmmObj.NComponents, ...
							'NDimensions', gmmObj.NDimensions, ...
							'NlogL', []);

		[~,gmmObjFake.NlogL] = cluster(gmmObj, data);

		mdl =  fitnessFunc('mdl',gmmObjFake, nObj);
	else 
		if isfield(resStruct,'fitness')
			mdl = resStruct.fitness;
		else 
			if isfield(resStruct,'m')
				mdl =  fitnessFunc('mdl', resStruct.m, nObj);
			else
				mdl =  fitnessFunc('mdl', resStruct, nObj);
			end
		end
	end
end
