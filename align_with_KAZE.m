function tform=align_with_KAZE(I1,I2)
I1=I1./nanmax(I1,[],'all');
I2=I2./nanmax(I2,[],'all');

% points1 = detectSURFFeatures(I1,'MetricThreshold',200);
% points2 = detectSURFFeatures(I2,'MetricThreshold',200);
points1 = detectKAZEFeatures(I1);
points2 = detectKAZEFeatures(I2);

[features1,valid_points1] = extractFeatures(I1,points1);
[features2,valid_points2] = extractFeatures(I2,points2);

indexPairs = matchFeatures(features1,features2,'MatchThreshold',50.000000,'MaxRatio',0.500000);

if isempty(indexPairs)
    tform=[];
else
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
tform=estimateGeometricTransform2D(matchedPoints2,matchedPoints1,'affine');
end