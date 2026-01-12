info = imfinfo(filename);
numFrames = numel(info);

imgStack = zeros(info(1).Height, info(1).Width, numFrames, 'uint16');

for k = 1:numFrames
    imgStack(:,:,k) = double(imread(filename, k));
    meanVal(k,1) = mean(double(imgStack(:,:,k)), 'all');
    intStack(:,k)=reshape(imgStack(:,:,k),[],1);
end

medVal = squeeze(median(median(imgStack,1),2));

interactiveImageStack(imgStack)