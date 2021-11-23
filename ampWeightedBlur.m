function out = ampWeightedBlur(in,weight,sigma)
% performs an amplitude weighted blur (guassian blur)
% weight is typically the amplitude obtained from a steerable complex
% filter, and in is phase

% weight = weight + eps; % apparently a trick to reduce noise claimed by the authors
% out = imgaussfilt(in.*weight,sigma);
% weightMat = imgaussfilt(weight,sigma);

out = zeros(size(in));
weightMat = zeros(size(in));
w = [1 6 15 20 15 6 1]/64;
for i = 1:size(in,3)
    out(:,:,i) = conv2(w,w,in(:,:,i).*weight(:,:,i),'same');
    weightMat(:,:,i) = conv2(w,w,weight(:,:,i),'same');
end

%%% alternative implementation:
% kernel = fspecial('gaussian',ceil(4*sigma),sigma);
% out = imfilter(in.*weight,kernel,'circular');
% weightMat = imfilter(weight,kernel,'circular');
%%%

out = out./weightMat;

end

