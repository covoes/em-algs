function nlogl = getNLOGL(resStruct, r, t, heldout, data)
	if isempty(resStruct) || isempty(resStruct{r})
		nlogl = NaN;
		return
	end
	if isempty(resStruct{r}{t})
		if t > 1 
			for t2=t-1:-1:1 
				if ~isempty(resStruct{r}{t2})
					nlogl = getNLOGLIntern(resStruct{r}{t2},heldout, data);
					return
				end
			end
			nlogl = NaN;
		else
			nlogl = NaN;
		end
	else
		nlogl = getNLOGLIntern(resStruct{r}{t},heldout, data);
	end
end

function nlogl = getNLOGLIntern(resStruct, heldout, data)
	if isa(resStruct, 'gmdistribution')
		%MREM
		[~,nlogl] = posterior(resStruct, heldout);
	else
		if isfield(resStruct,'determinant')
			%ESMEM
			[~,nlogl] = posterior(resultsToGMM(resStruct), heldout);
		else 
			if isfield(resStruct,'bits')
				%GAEM
				[~,nlogl] = posterior(resultsGAEMToGMM(resStruct), heldout);
			else 
				if isfield(resStruct, 'medoids')
					%GAKREM
					[~,nlogl] = posterior(resultsGAKREMToGMM(resStruct, data), heldout);

				else
					if isfield(resStruct,'m')
						%OMREM com contador de loop
						[~,nlogl] = posterior(resStruct.m, heldout);
					end
				end
			end
		end
	end
end
