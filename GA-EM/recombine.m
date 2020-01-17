function [ childrens ] = recombine( P, data, maxClusters, pRecombination )
%RECOMBINE Perform single-point crossover in two random individuals to
%generate H childrens (H=pRecombination*populationSize)

    populationSize = size(P, 1);
    %define number of childrens
    H = round( pRecombination * populationSize );
    childrens = initialize( H, data, maxClusters );

	for i=1:2:H


		%select parents randomly and point of crossover
		parents = randsample( populationSize, 2 );
		point = randi(maxClusters);

		%copy data to childrens
		%children gets the first part from the father and the second from the mother        
		childrens(i).bits = [ P( parents(1) ).bits( 1:point-1 ) P( parents(2) ).bits( point:maxClusters ) ];

		%if i==H the last children will be discarded
		if( i<H )
			%other children gets the first part from the mother and the second from the father        
			childrens(i+1).bits = [ P( parents(2) ).bits( 1:point-1 ) P( parents(1) ).bits( point:maxClusters ) ];
		end

		for k=1:maxClusters
			%using 0,1 to identify parent to use complement after
			if( k < point)
				parentToCopy = 0;
			else 
				parentToCopy = 1;
			end
			childrens(i).clusters(k).mean = P( parents( parentToCopy+1 ) ).clusters(k).mean;
			childrens(i).clusters(k).covariance = P( parents( parentToCopy+1 ) ).clusters(k).covariance;

			%if i==H the last children will be discarded
			if( i<H )
				%the second children must get the info from the other father that the first did (thats why the ~)
				childrens(i+1).clusters(k).mean = P( parents( ~parentToCopy+1 ) ).clusters(k).mean;
				childrens(i+1).clusters(k).covariance = P( parents( ~parentToCopy+1 ) ).clusters(k).covariance;
			end
		end
			
	end
end

