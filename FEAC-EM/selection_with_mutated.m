function [P] = selection_with_mutated(P,Pmut,numElitists);
%SELECTION Performs elitism and roulette wheel selection on the population formed by parents and mutated individuals
%
%Assumes that the intent is to minimize the fitness function

global DEBUG;

	fullPop = [ P'; Pmut' ]; 
	fitness = [fullPop(:).fitness];
	sizePopulation = length(P);

	%elitism
	[bestFit idxSorted] = sort(fitness,'ascend');
	best = idxSorted(1:numElitists);
	%roulette
	ranks = zeros(length(fullPop),1);
	ranks(idxSorted) = length(fullPop):-1:1;
	chance = ranks/sum(ranks);	
	chosen = roulette( chance, sizePopulation-numElitists );
	if DEBUG
		fprintf(DEBUG,'#SELECTION\nFITNESS:\n\t%s',mat2str(fitness,4));
		fprintf(DEBUG,'\nNUM K:\n\t%s',mat2str([P(:).numClusters],4));
		fprintf(DEBUG,'\nPROBS:\n\t%s',mat2str(chance,4));
		fprintf(DEBUG,'\nELITISM:\n\t%s',mat2str(best));
		fprintf(DEBUG,'\nCHOSEN:\n\t%s\n',mat2str(chosen,4));
	end

	P(1:numElitists) = fullPop(best);
	P(numElitists+1:sizePopulation) = fullPop( chosen );

end
