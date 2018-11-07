function [neurons] = SOM(dataset,mapWidth,maxIteration,startLearningRate)
[row col]=size(dataset);
Nofeatures = col; 
dataset = dataset(1:row,1:Nofeatures);
neurons = rand(mapWidth,Nofeatures);
%neurons = dataset(1:mapWidth,1:Nofeatures);
mapRadius = max(mapWidth, 1) / 2;
timeConstant = maxIteration / log(mapRadius);
i = 1;
group = zeros(row,1);
for i = 1:maxIteration
    randval = randi(row-1) + 1;
    vector = dataset(randval,1:Nofeatures);
    newmatrix = [neurons;vector];
    distance = pdist(newmatrix);
    s = squareform(distance);
    [mindist, bmu] = getMinimum(s, size(s,1));
    group(randval,1) = bmu;
    neighbourhoodRadius =  mapRadius * exp ( -i / timeConstant);
    radiusSquare = neighbourhoodRadius * neighbourhoodRadius;
    learningRate = startLearningRate * exp (-i / maxIteration);
 for n = 1:size(neurons,1)
        temp = pdist([getGridPosition(n,mapWidth); getGridPosition(bmu,mapWidth)]);
        distSquare = temp(1,1);
      if(distSquare < radiusSquare)
             m_dInfluence = exp(-(distSquare) / (2*radiusSquare));
             for c = 1:length(neurons(n,:))
                 neurons(n,c) = neurons(n,c) + m_dInfluence * learningRate * (vector(c) - neurons(n,c));
             end
      end
 end
 
end



end
 



 

 
