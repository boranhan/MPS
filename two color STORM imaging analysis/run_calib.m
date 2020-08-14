clear all
close all

% calibrate in 3D the two channels.
% generate the warping matrix from Cy5 channel to the Cy3 channel.

FileName = ['calib/movie_0004'];
[ImageStack, InfoFile] = ReadZStack(FileName,111,5);
for i = 1:size(ImageStack, 3)
    LeftImage(:,:,i) = ImageStack(:, 1:256, i);
    RightImage(:,:,i) = ImageStack(:, 257:512, i);
    StdLeft(i) = std2(LeftImage(:,:,i));
    StdRight(i) = std2(RightImage(:,:,i));
end

ImageMeanLeft = mean(LeftImage,3);
ImageMeanRight = mean(RightImage,3);
% figure(100)
% subplot(1,2,1)
% imagesc(ImageMeanLeft)
% colormap gray
% axis equal
% subplot(1,2,2)
% imagesc(ImageMeanRight)
% colormap gray
% axis equal

RGB(:,:,1) = ImageMeanLeft/max(max(ImageMeanLeft))*2;
RGB(:,:,2) = ImageMeanRight/max(max(ImageMeanRight))*2;
RGB(:,:,3) = zeros(256,256);
RGB(find(RGB>1)) = 1;
figure(200)
imagesc(RGB)
axis equal

WhetherROI = questdlg('Do you want to select ROIs ?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
    selectROI;    
    lmol_match =[];
    rmol_match =[];
    for i = 1:length(roiList)
        [Xfit, Yfit, Zfit] = fitFoci(LeftImage, roiList(i),2*i-1,1);      
        lmol_match = cat(1, lmol_match, [Xfit, Yfit, Zfit]);
        [Xfit, Yfit, Zfit] = fitFoci(RightImage, roiList(i),2*i,1);   
        rmol_match = cat(1, rmol_match, [Xfit, Yfit, Zfit]);
    end
    tform = cp2tform(lmol_match(:,1:2),rmol_match(:,1:2),'projective');
    save('tform.mat','tform');
end

GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
StartPoint = [max(StdLeft) 10 1 0];
f_left = fit([1:length(StdLeft)]', StdLeft', GaussEqu, 'Start', StartPoint);
StartPoint = [max(StdRight) 10 1 0];
f_right = fit([1:length(StdRight)]', StdRight', GaussEqu, 'Start', StartPoint);
figure(300)
subplot(1,2,1)
plot(f_left, 1:length(StdLeft), StdLeft);
legend('off');
subplot(1,2,2)
plot(f_right, 1:length(StdRight), StdRight);
legend('off');
DeltaZ = f_left.b-f_right.b;
save('DeltaZ.mat','DeltaZ');

TransImg = imtransform(RGB(:,:,1), tform, 'XData', [1 256], 'Ydata', [1 256]);
RGB(:,:,1) = TransImg;
figure(500)
imagesc(RGB)
axis equal

