close all
clear all
fclose all;
ipath = '';
filename = dir(fullfile(strcat([ipath '\'],'*g.spe')));
min_dist_total=[];
ratio=[];
for ii=1:length(filename)
    Im_green = readSPE([ipath '\'], filename(ii).name);
    
    Im_green=mean(Im_green,3);
    Im_green=Im_green-min(Im_green(:));
%     Im_green = imresize(Im_green,5);
    Im_red= readSPE([ipath '\'], [filename(ii).name(1:end-5) 'r.spe']);
    Im_red=mean(Im_red,3);
    Im_red=Im_red-min(Im_red(:));
%     Im_red = imresize(Im_red,5);
    backgroundIm_1 = imopen(Im_red, strel('disk', 20));
    Im_red_contrast = Im_red - backgroundIm_1;
    
    Im_red_contrast = wiener2(Im_red_contrast,[3 3]);
    Im_red_contrast= imadjust(Im_red_contrast/max(Im_red_contrast(:)));
    Im_green= imadjust(Im_green/max(Im_green(:)));
%     backgroundIm_1 = imopen(Im_red_contrast, strel('disk', 10));
%     Im_red_contrast = Im_red_contrast - backgroundIm_1;
%     backgroundIm_2 = imopen(Im_green, strel('disk', 20));
%     Im_green = Im_green - backgroundIm_2;
    
    figure(5);imshow(Im_red_contrast)
    image_binary_green = double(imbinarize(Im_green,graythresh(Im_green)*0.7));
%     image_binary_green = double(imbinarize(Im_green/max(Im_green(:)),'adaptive','Sensitivity',0.01));
    image_binary_green = imdilate(image_binary_green, strel('disk', 1));
    image_binary_red = double(imbinarize(Im_red_contrast/max(Im_red_contrast(:)),'adaptive','Sensitivity',0.65));%graythresh(Im_red_contrast/max(Im_red_contrast(:)))*0.8));
    image_binary_green= double(bwareaopen(image_binary_green, 1000,4));
%     image_binary_red= double(bwareaopen(image_binary_red, 100,4));
    renderedStack=cat(3,image_binary_red,image_binary_green);
    renderedStack=cat(3,renderedStack,zeros(size(renderedStack,1),size(renderedStack,2)));
%     figure(4);imshow(renderedStack, []);
    image_binary_red=double(image_binary_red)-double(image_binary_green);
    image_binary_red(image_binary_red<0)=0;
    image_binary_red= double(bwareaopen(image_binary_red, 500,4));
    
%     image_binary_area = bwareaopen(image_binary_area, 2);
    figure(2);imshow(Im_green/max(Im_green(:)));
    figure(1);imshow(image_binary_red);
    image_binary_red_distance = bwdist(image_binary_red);
%     figure(3);imshow(image_binary_area);
    image_intensity_red_lbl=bwlabel(image_binary_red);
    min_dist_total(ii)=mean(image_binary_red_distance(image_binary_green==1));
%     for iii=1:max(image_intensity_red_lbl(:))
%         index_red=find(image_intensity_red_lbl==iii);
%         min_dist=min(image_binary_green_distance(index_red));
%         min_dist_total=[min_dist,min_dist_total];
%     end
    image_binary_green_dilate = imdilate(image_binary_green, strel('disk', 3))-image_binary_green;
    renderedStack=cat(3,image_binary_red,image_binary_green_dilate);
    renderedStack=cat(3,renderedStack,zeros(size(renderedStack,1),size(renderedStack,2)));
    figure(4);imshow(renderedStack, []);
    ratio(ii)=sum(image_binary_red(image_binary_red&image_binary_green_dilate))/sum(image_binary_green_dilate(:));%/sum(image_binary_red(:));%/sum(image_binary_green(image_binary_red|image_binary_green));
%     ratio(ii)=sum(image_binary_red(:));%/sum(image_binary_green(image_binary_red|image_binary_green));
end
