function constraintList = constraintListOfMatrix(constraints,numObj)
%CONSTRAINTLISTOFMATRIX Generates a list of each constraint an object is related to. Each entry on the list is the key for the matrix 


objs = constraints(:,1:2);
objs = unique(objs(:));
constraintList = cell(1, numObj);
for o=objs'
	constraintList{o}=find( (constraints(:,1) == o) | (constraints(:,2) == o) );
end
