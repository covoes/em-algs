function [P] = mutation(P, numElitists, splitFName, splitOpName, mergeFName, mergeOpName, maxKMSIter, maxClusters, data)
%MUTATION Performs mutation on the indivuals in the population, except the ones mantained by elitism
%
%It assumes that the first individuals are the ones mantained by elitism
%Assumes that is minimizing the fitness function

global DEBUG;


for i=numElitists+1:length(P)

	mutOpToApply = NaN;
	%verify that the current operation (split|merge) is improving the individual
	if (P(i).fitness - P(i).lastFitness < 0 )
		mutOpToApply = P(i).lastMutOperator;
		%no mutation was performed last generation in this individual
		if( isnan(mutOpToApply) )
			mutOpToApply = randi([0 1]);
		end
	else
		%change the current mutation operation (split|merge) 
		%if nothing was done, select randomly
		if( isnan(P(i).lastMutOperator) )
			mutOpToApply = randi([0 1]);
		else
			mutOpToApply = ~P(i).lastMutOperator;
		end
	end


	if DEBUG
		MUTNAME='merge';
		if mutOpToApply == 0
			MUTNAME='split';
		end
				
		fprintf(DEBUG,'#MUTATION\nIND(%d)\nLAST OPERATOR: %d\nLAST FITNESS %.4f\nOPTOAPPLY %s\nCURRENT FITNESSS %.4f\n', i,P(i).lastMutOperator, P(i).lastFitness, MUTNAME, P(i).fitness);
		fprintf(DEBUG,'#MUTATION\nOLD INDIVIDUAL (%d):%s\n',i,info_individual(P(i)));
	end


	if mutOpToApply == 0
		P(i)=split(P(i), splitFName, splitOpName, maxKMSIter, maxClusters, data);
	else
		P(i)=merge(P(i), mergeFName, mergeOpName, data);
	end

	%update individual
	P(i).lastMutOperator = mutOpToApply;

	if DEBUG
		fprintf(DEBUG,'#MUTATION\nNEW INDIVIDUAL (%d):%s\n',i,info_individual(P(i)));
	end


end
