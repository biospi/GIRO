function dispOverlay(img1, img2, RT, MZ)

tmp = zeros(size(img1,1), size(img1,2), 3);

% tmp(:,:,1) = log(img1+10) ./ max(max(log(img1+10)));
% tmp(:,:,2) = log(img2+10) ./ max(max(log(img2+10)));
% tmp(:,:,3) = log(img1+10) ./ max(max(log(img1+10)));

tmp(:,:,1) = (img1+1) ./ max(max((img1+1)));
tmp(:,:,2) = (img2+1) ./ max(max((img2+1)));
tmp(:,:,3) = (img1+1) ./ max(max((img1+1)));

tmp = 1 - tmp;

figure;
imagesc(MZ,RT,tmp);
xlabel('M/Z, Thompson')
ylabel('Retention Time, min')
title('LC-MS')
