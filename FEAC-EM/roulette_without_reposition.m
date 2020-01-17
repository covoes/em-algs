function [chosen] = roulette_without_reposition( probs, numToChoose )
%ROULETTE_WITHOUT_REPOSITION Performs roulette wheel selection without replacement

chosen = zeros(numToChoose,1);	

for i = 1:numToChoose
	selected = roulette( probs, 1 );
	probs = probs / (1-probs(selected));
	probs( selected ) = 0;
	chosen(i) = selected;
end

end
