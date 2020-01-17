function [X,W,M,Co,lbls,sepm] = csepGenerator(n,k,d,c,e,unif,varsep)
%mixgen - Gaussian mixture generator 
%
%[X,T,L1,L2] = mixgen(n,m,k,d,c,e) 
%  n - size of training set 
%  m - size of test set 
%  k - number of components
%  d - dimension
%  c - separation degree
%  e - maximum eccentricity
%  varsep - bound for variation i.e. varsep = 0.1 mean that each pair of gaussian can be c-.1 or c+.1 apart.
%returns
%  X  - training set (n x d)
% W...weight matrix
% M...Mean matrix
% Co...Covarianz matrix
% lbls...Label vector
% Nikos Vlassis, 2000
% modified by pernkopf
% for definitions see (Dasgupta, 1999)
% modified by covoes, separation now has a bound

	R=zeros(k,d^2);

	% mixing weights
	if (unif==1)
		W=ones(k,1)*1/k;
	else
		while 1
			W = rand(k,1); 
			W = W / sum(W);
			if all(W > 1/(2*k)) 
				break;
			end
		end
	end

	% create c-separated Gaussian clusters of maximum eccentricity e
	trials = 1;

	while 1
		count=0;    
		X = [];
		M = randn(k,d)*sqrt(k)*sqrt(c)*trials;
		Trace = zeros(k,1);
		Me=zeros(k,d);
		Co=zeros(d,d,k);
		lbls=zeros(n,1);
		lstIdx=0;

		for j = 1:k
			count=count+1;  
			U = rand(d,d)-0.5; 
			U = sqrtm(inv(U*U')) * U;
			L = diag(rand(d,1)*(e-1)+1).^2;
			msg = 1;
			while msg
				[C,msg] = chol(U*L*U');
			end
			R(j,:)=C(:)';

			nj = ceil(n*W(j));
			Xj = randn(nj,d) * C;
			X = [X; repmat(M(j,:),nj,1) + Xj];
			Trace(j) = trace(cov(Xj));
			Me(count,:)=mean(Xj)'+M(j,:)';
			Co(:,:,count)=cov(Xj);
			W(count)=size(Xj,1)/n;
			lbls(lstIdx+1:lstIdx+size(Xj,1)) = count;
			lstIdx = lstIdx+size(Xj,1);

		end
		

		for r=1:10 
			% try to fix degree of separation once, and then select the minimum csep, to be equal to the wanted
			for i = 1:k-1
				for j = i+1:k
					%if norm(M(i,:)-M(j,:)) < c * sqrt(2*max(max(eig(Co(:,:,i))),max(eig(Co(:,:,j)))))
					%if norm(M(i,:)-M(j,:)) < c * sqrt(max(trace(Co(:,:,i)),trace(Co(:,:,j))))
					%if norm(M(i,:)-M(j,:)) < c * sqrt(max(Trace(i),Trace(j)))
					%new part, that adjust c-separation
					while true
						sep = norm(M(i,:) - M(j,:)) / (sqrt(d * max([max(eig(Co(:,:,i))) max(eig(Co(:,:,j))) ]))); 
						alpha = 0.01;
						if ( sep > (c+varsep) )
							%fprintf('fixing (approaching) %d : %d -> %g\n',i,j,sep);
							%M(i,:) = M(i,:) + alpha *(M(j,:) - M(i,:));
							M(j,:) = M(j,:) + alpha *(M(i,:) - M(j,:));
						elseif( sep < c )
							%fprintf('fixing (moving away) %d : %d -> %g\n',i,j,sep);
							%M(i,:) = M(i,:) - alpha *(M(j,:) - M(i,:));
							M(j,:) = M(j,:) - alpha *(M(i,:) - M(j,:));
						else
							%fprintf('fixed %d : %d -> %g\n',i,j,sep);
							break;
						end
					end
				end
			end
			sepm = nan(k);
			for i = 1:k-1
				for j = i+1:k
					sepm(i,j) = norm(M(i,:) - M(j,:)) / (sqrt(d * max([max(eig(Co(:,:,i))) max(eig(Co(:,:,j))) ]))); 
				end
			end

			[i j] = find(sepm == min(sepm(:)));

			while true
				sep = norm(M(i,:) - M(j,:)) / (sqrt(d * max([max(eig(Co(:,:,i))) max(eig(Co(:,:,j))) ]))); 
				alpha = 0.01;
				if ( sep > (c+0.005) )
					%fprintf('fixing (approaching) %d : %d -> %g\n',i,j,sep);
					%M(i,:) = M(i,:) + alpha *(M(j,:) - M(i,:));
					M(j,:) = M(j,:) + alpha *(M(i,:) - M(j,:));
				elseif( sep < c )
					%fprintf('fixing (moving away) %d : %d -> %g\n',i,j,sep);
					%M(i,:) = M(i,:) - alpha *(M(j,:) - M(i,:));
					M(j,:) = M(j,:) - (alpha/10) *(M(i,:) - M(j,:));
				else
					%fprintf('fixed %d : %d -> %g\n',i,j,sep);
					break;
				end
			end

		end

		%if fixing one, create a problem with the other.. try the whole thing again
		errorr = 0;
		for i = 1:k-1
			for j = i+1:k
				sep = norm(M(i,:) - M(j,:)) / (sqrt(d * max([max(eig(Co(:,:,i))) max(eig(Co(:,:,j))) ]))); 
				if( abs(sep-c) > varsep )
					%fprintf('PROBLEM %d : %d -> %g\n',i,j,sep);
					errorr = 1; 
					break;
				end
			end
		end

		if ~errorr
			break;
		end
		trials = trials + 1;
	end


	for i = 1:k-1
		for j = i+1:k
			sepm(i,j) = norm(M(i,:) - M(j,:)) / (sqrt(d * max([max(eig(Co(:,:,i))) max(eig(Co(:,:,j))) ])));
		end
	end
	sepm

	[X lbls] = random( gmdistribution( M, Co, W), n);
end
