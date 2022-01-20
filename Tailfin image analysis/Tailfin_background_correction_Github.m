% 
% ************************************************************************
% (c) 2019, Nathaniel R. Campbell, Xavier and White Labs, MSKCC
% campben2@mskcc.org | nathaniel.r.campbell@gmail.com
% 
% make binary image from ratio of GFP and tdTomato
% ************************************************************************

% Requires:
% - bfopen.m (bioformats plugin)
%

%% open files

fileName = '/Users/weiss/Desktop/SOS MEK tailfin imaging 051221/Refametinib + SOSi/Snap-00083.czi';

Im = bfopen(fileName);


%% show images


figure
subplot(1,3,1)
imagesc(Im{1,1}{1,1})
title('BF')


subplot(1,3,2)
imagesc(Im{1,1}{4,1})
colorbar
title('tdT')

subplot(1,3,3)
imagesc(Im{1,1}{7,1})
colorbar
title('GFP')


%% gaussian background correction
% information purposes only

Ig = imgaussfilt(Im{1,1}{7,1},100);


figure
subplot(1,3,1)
imagesc(Im{1,1}{6,1})
title('original')

subplot(1,3,2)
imagesc(Ig)
title('gaussian blur')

subplot(1,3,3)
imagesc(imsubtract(Im{1,1}{7,1}, Ig));
title('gaussian corrected')




%% subtract median
% information purposes only

Ig = imgaussfilt(Im{1,1}{7,1},100);


figure
subplot(1,2,1)
imagesc(Im{1,1}{7,1})
colorbar
title('original')


ImMed = median(Im{1,1}{7,1}(:));


subplot(1,2,2)
imagesc(Im{1,1}{7,1} - double(ImMed));
title('median corrected')


%% histogram of GFP
% information purposes only

figure
histogram(double(Im{1,1}{7,1}(:)))
xlabel('GFP Intensity')
ylabel('Count')


%% ratio with gaussian correction
% information purposes only

ImDG = imgaussfilt(double(Im{1,1}{7,1}) ./ double(Im{1,1}{4,1}),100);


figure
subplot(1,3,1)
imagesc(double(Im{1,1}{7,1}) ./ double(Im{1,1}{4,1}))
colorbar
title('GFP/tdT ratio')

subplot(1,3,2)
imagesc(ImDG)
colorbar
title('GFP/tdT ratio gaussian')

subplot(1,3,3)
imagesc(imsubtract(double(Im{1,1}{7,1}) ./ double(Im{1,1}{4,1}), ImDG))
colorbar
title('GFP/tdT ratio gaussian corrected')


%% ratio median (use me!!)

ratThresh = 800; % threshold for segmentation, adjust as necessary, 3000 was original




ImRat = imresize(1000.*(double(Im{1,1}{7,1}) ./ double(Im{1,1}{4,1})), 0.5);
ImDG = imgaussfilt(ImRat,100);



figure
subplot(1,3,1)
imagesc(ImRat)
colorbar
title('ratio image')

subplot(1,3,2)
imagesc(ImRat - double(median(ImRat(:))))
colorbar
title('ratio median subtract')

subplot(1,3,3)
imagesc(im2bw(uint16(ImRat - double(median(ImRat(:)))), ratThresh/(2^16-1)))
title('binarize')


%% histogram of ratio median corrected (zoomed in on y-axis)

ImRatCor = ImRat - double(median(ImRat(:)));

figure
histogram(double(ImRatCor(:)), 'normalization', 'probability')
xlabel('GFP/tdT corrected')
ylabel('Frequency')

set(gca, 'ylim', [0, 0.005])

%% save data as tif stack

ImBinary = imresize(im2bw(uint16(ImRat - double(median(ImRat(:)))), ratThresh/(2^16-1)),2);

writeLoc = '~/Desktop'; % where do you want to write the file



imwrite(Im{1,1}{1,1}, [writeLoc '/' fileName(end-13:end-4) '_binary.tif'], 'compression', 'none')
imwrite(uint16(imresize(ImRat,2)), [writeLoc '/' fileName(end-13:end-4) '_binary.tif'], 'compression', 'none', 'WriteMode', 'append')
imwrite(uint16(ImBinary), [writeLoc '/' fileName(end-13:end-4) '_binary.tif'], 'compression', 'none', 'WriteMode', 'append')



