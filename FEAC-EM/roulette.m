function [chosen] = roulette( probs, numToChoose )
%ROULETTE Performs roulette wheel selection with replacement

chosen = zeros(numToChoose,1);	

for i = 1:numToChoose
	draw = rand;
	acc = probs(1);
	chosen(i) = 1;
	while (draw > acc)
		chosen(i) = chosen(i) + 1;
		acc = acc + probs(chosen(i)); 
	end
end

end
