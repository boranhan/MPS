ipath = '';
filename = dir(fullfile(strcat([ipath '\'],'*g.spe')));
cell_distance=[];
scale=5;
step_size=10;
crop_size=30;
for ii=1:15%length(filename)
    distance=[];
    Im_green = readSPE([ipath '\'], filename(ii).name);
    Im_green=mean(Im_green,3);
    Im_green=Im_green-min(Im_green(:));
%     Im_green = imopen(Im_green, strel('disk', 1));
    Im_red= readSPE([ipath '\'], [filename(ii).name(1:end-5) 'r.spe']);
    Im_red=mean(Im_red,3);
    Im_red=Im_red-min(Im_red(:));
    Im_red=imresize(Im_red,scale);
    Im_green=imresize(Im_green,scale);
    
    image_binary_green = double(imbinarize(Im_green/max(Im_green(:)),graythresh(Im_green/max(Im_green(:)))*1.2));
    image_binary_green = imdilate(image_binary_green, strel('disk', 40));
    image_binary_red = double(imbinarize(Im_red/max(Im_red(:)),graythresh(Im_red/max(Im_red(:)))*1));
    image_binary_green= double(bwareaopen(image_binary_green, 4000*scale,4));
    image_binary_red= double(bwareaopen(image_binary_red, 100*scale,4));
    renderedStack=cat(3,image_binary_red,image_binary_green);
    renderedStack=cat(3,renderedStack,zeros(size(renderedStack,1),size(renderedStack,2)));
    figure(4);imshow(renderedStack, []);

    image_binary_red=double(image_binary_red)-double(image_binary_green);
    image_binary_red(image_binary_red<0)=0;
    image_binary_red= double(bwareaopen(image_binary_red, 5000*scale,4));
    
%     image_binary_area = bwareaopen(image_binary_area, 2);
    figure(2);imshow(image_binary_green);
    figure(5);imshow(image_binary_red);
    image_binary_red_boundary=zeros(size(image_binary_green,1),size(image_binary_green,2));
    B = bwboundaries(image_binary_red);
    for k = 1:length(B)
        boundary = B{k};
        if length(boundary)>3
            for kk = 1:length(boundary)
                image_binary_red_boundary(boundary(kk,1), boundary(kk,2))=1;
            end
        end
    end
    
    imshow(image_binary_red_boundary)
    image_binary_red_boundary = logical(image_binary_red_boundary);
    image_binary_red_boundary_lbl=bwlabel(image_binary_red_boundary);
%     imagesc(image_binary_red_boundary_lbl);
    for k = 100:step_size*scale:size(image_binary_green,1)-2*crop_size*scale
        for kk = 100:step_size*scale:size(image_binary_green,2)-2*crop_size*scale
            image_crop_boundary=image_binary_red_boundary(k:k+crop_size*scale,kk:kk+crop_size*scale);
            image_crop_boundary = logical(image_crop_boundary);
            image_crop_boundary= double(bwareaopen(image_crop_boundary, 5*scale));
            image_crop_boundary_lbl=bwlabel(image_crop_boundary);
            image_crop=image_binary_red(k:k+crop_size*scale,kk:kk+crop_size*scale);
            image_crop = imdilate(image_crop, strel('disk', 2));
            image_crop_lbl=bwlabel(image_crop);
            
%             figure(3);imagesc(image_crop_boundary_lbl);
%             figure(4);imagesc(image_crop_lbl);
            distance_region=[];
            delete=0;
            
            for kkkk = 1:max(image_crop_lbl(:))
                boundary_number=unique(image_crop_boundary_lbl(image_crop_lbl==kkkk));
                if length(boundary_number)>1
                    for kkk = boundary_number(2:end)'
                        image_crop_1=zeros((crop_size*scale+1)*(crop_size*scale+1),1);
                        index=find(image_crop_boundary_lbl==kkk);
                        [col,row]=ind2sub(crop_size*scale+1,index);
                        if ~((max(col)==(crop_size*scale+1)||min(col)==1)&&(max(row)==(crop_size*scale+1)||min(row)==1))&&~(max(col)==(crop_size*scale+1)&&min(col)==1)&&~(max(row)==(crop_size*scale+1)&&min(row)==1)
                            delete=0;
                        end
                        image_crop_1(index)=1;
                        image_crop_1=reshape(image_crop_1,[],(crop_size*scale+1));
                        D = bwdist(image_crop_1);   
                        for kkk_1 = boundary_number(2:end)'
                            distance_region=[distance_region,mink(D(image_crop_boundary_lbl==kkk_1),1)];
                        end
                    end
                end
            end

            if delete==1
                distance_region=[];
            end
            [C,idx]=unique(distance_region,'stable');%don't sort
            idx=setxor(idx,1:numel(distance_region));
            repeating_values=distance_region(idx);
            distance=[distance,unique(distance_region(distance_region>0))];
        end
    end
    cell_distance=[cell_distance,mean(distance(distance>0)/scale)];
end

min_dist_total=cell_distance(~isnan(cell_distance));

