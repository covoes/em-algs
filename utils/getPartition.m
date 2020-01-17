function idx = getPartition(resStruct, r, t, data)
	nObj = size(data,1);
	idx = ones(1,nObj);
	if r == -1
		gmm = getObj(resStruct);
	elseif isempty(resStruct) || isempty(resStruct{r})
		return;
	else
		gmm = [];
		if isempty(resStruct{r}{t})
			if t > 1
				%find the last one 
				for t2=t-1:-1:1
					if	~isempty(resStruct{r}{t2})
						gmm = getObj(resStruct{r}{t2});
					end
				end
			end
		else
			gmm = getObj(resStruct{r}{t});
		end
	end
	if isempty(gmm)
		return;
	end
	idx = cluster(gmm, data);
end


function gmm = getObj(resStruct)
	if isa(resStruct, 'gmdistribution')
		%MREM
		gmm = resStruct;
	else
		if isfield(resStruct,'determinant')
			%ESMEM
			gmm = resultsToGMM(resStruct);
		else 
			if isfield(resStruct,'bits')
				%GAEM
				gmm = resultsGAEMToGMM(resStruct);
			else isfield(resStruct,'m')
				gmm = resStruct.m;
			end
		end
	end
end

