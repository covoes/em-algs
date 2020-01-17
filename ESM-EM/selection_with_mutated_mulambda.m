function [P] = selection_with_mutated_mulambda(P,Pmut,sizePopulation);
%SELECTION Performs(mu+lambda) selection 
%
%Assumes that the intent is to minimize the fitness function

global DEBUG;

	fullPop = [ P'; Pmut' ]; 
	fitness = [fullPop(:).fitness];

	%elitism
	[bestFit idxSorted] = sort(fitness,'ascend');
	P = fullPop(idxSorted(1:sizePopulation));

end
