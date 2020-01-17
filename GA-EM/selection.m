function [ P ] = selection( parents, childrens, populationSize)
%SELECTION Implement (K,H)-strategy, which selects populationSize from the
%population formed by the parents and childrens population
%
%Parameters:
%parents       : parents population
%childrens     : childrens population
%populationSize: size of the next population

    P = [ parents; childrens ];
    
    [MDLs, idx] = sort( [ P.fitness ], 2, 'ascend' );

    P = P(idx(1:populationSize));

end

