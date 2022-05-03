clear all
close all

% calibrate in 3D the two channels.
% generate the warping matrix from Cy5 channel to the Cy3 channel.
dataPath = '';
FileName = ['colocalize_0011_5'];
[ImageStack, InfoFile] = ReadZStack(FileName,10,1);
append=zeros(512,512,9);
for i = 1:size(ImageStack, 3)
ImageStack(:,:,i) = rot90(flipud(ImageStack(:,:,i)));
end
ImageStack=cat(1,ImageStack,append);
for i = 3:size(ImageStack, 3)
    LeftImage(:,:,i-2) = ImageStack(1:256,253:508, i);
    RightImage(:,:,i-2) = ImageStack(254:509, 252:507, i);
    StdLeft(i-2) = std2(LeftImage(:,:,i-2));
    StdRight(i-2) = std2(RightImage(:,:,i-2));
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

RGB(:,:,1) = ImageMeanLeft/max(max(ImageMeanLeft))*1.2;
RGB(:,:,2) = ImageMeanRight/max(max(ImageMeanRight))*1.5;
RGB(:,:,3) = zeros(256,256);
RGB(find(RGB>1)) = 1;
figure(200)
chim=[0.99 1];
imagesc(RGB,chim)
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
%% 

[Datafiles, InfoFile] = ReadDax([FileName, '.dax'],'startFrame', 1, 'endFrame', 510);
Datafiles=double(Datafiles);
% for i = 1:size(Datafiles, 3)
% Datafiles(:,:,i) = rot90(flipud(Datafiles(:,:,i)));
% end
append=zeros(512,512,9);
ImageStack=cat(1,ImageStack,append);

for i = 1:size(Datafiles, 3)
    LeftDatafiles(:,:,i) = Datafiles(1:256,253:508, i);
    RightDatafiles(:,:,i) = Datafiles(254:509, 252:507, i);
    TransImg = imtransform(LeftDatafiles(:,:,i), tform, 'XData', [1 256], 'Ydata', [1 256]);
    LeftDatafiles(:,:,i) = TransImg;
end
LeftDatafiles=uint16(LeftDatafiles);
InfoFile.vstart=1;
InfoFile.hstart=1;
InfoFile.frame_dimensions=[256 256];
newInfFile=InfoFile;
filteredDataPath = [dataPath 'colocalize_0011_5_left\'];
newInfFile.localPath = filteredDataPath;
WriteDAXFiles(LeftDatafiles, newInfFile);

RightDatafiles=uint16(RightDatafiles);
filteredDataPath = [dataPath 'colocalize_0011_5_right\'];
newInfFile.localPath = filteredDataPath;
WriteDAXFiles(RightDatafiles, newInfFile);
%% 
RightDatafiles=double(RightDatafiles)/max(max(max(double(RightDatafiles))))*2.7;
LeftDatafiles=double(LeftDatafiles)/max(max(max(double(LeftDatafiles))))*2.3;
for i = 1:size(Datafiles, 3)
    RGB_data(:,:,:,i) = cat(3,LeftDatafiles(:,:,i),RightDatafiles(:,:,i),zeros(256,256,1));
end
v = VideoWriter('colocalize_0011_5_high_contrast.avi');
open(v)
for i=1:size(RGB_data,4)
Img(:,:,:)=RGB_data(:,:,:,i);
Img(find(Img>1))=1;
figure(2);
Img=imadjust(Img,[0.1 1],[]);
imshow(Img,'InitialMagnification','fit');
mov(i)=getframe(gcf);
end
movie2avi(mov,'colocalize_0011_5_high_contrast.avi','compression', 'None');
%% 

for i=1:size(RGB_data,4)
Img_left(:,:,:)=RGB_data(:,:,:,i);
Img_left(:,:,2)=Img_left(:,:,1);
Img_left(:,:,3)=Img_left(:,:,1);
Img_left(find(Img_left>1))=1;
figure(2);
Img=imadjust(Img,[0.7 1],[]);
imshow(Img_left,'InitialMagnification','fit');
mov(i)=getframe(gcf);
end
movie2avi(mov,'colocalize_0011_5_high_contrast_left.avi','compression', 'None');

for i=1:size(RGB_data,4)
Img_right(:,:,:)=RGB_data(:,:,:,i);
Img_right(:,:,1)=Img_right(:,:,2);
Img_right(:,:,3)=Img_right(:,:,2);
Img_right(find(Img_right>1))=1;
Img=imadjust(Img,[0.3 1],[]);
figure(2);
imshow(Img_right,'InitialMagnification','fit');
mov(i)=getframe(gcf);
end
movie2avi(mov,'colocalize_0011_5_high_contrast_right.avi','compression', 'None');