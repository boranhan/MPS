%close all
ipath = 'R:\20201124 b2KD addKD WT dendritre comparison\msDIV21 addKD GFP-488 M2-Cy3 CTB-647 20x\';
filename = dir(fullfile(strcat([ipath '\'],'*g.spe')));
%distance = 0;
dist_array = [];
diameter_array = [];
w = 10;
ID_count = 0;
for ii=1:length(filename)
    Im_green = readSPE([ipath '\'], filename(ii).name);
    Im_green=mean(Im_green,3);
    Im_green=Im_green-min(Im_green(:));
    %imshow(Im_green,[])
    %hold on
    T_all = readtable([ipath filename(ii).name(1:end-4) '.txt'], 'ReadRowNames',true);
    Tracing_ID = unique(T_all.Tracing);
    for jj = 1: length(Tracing_ID)
        ID_count = ID_count+1;
        distance = 0;
        count = 0;
        T = T_all(strcmp(T_all.Tracing, Tracing_ID{jj}),:);%strcmp(t.a1, 'don'), :
        for j = 1:height(T)-1
            start_point = [T.X_pix_(j)+1 T.Y_pix_(j)+1];
            goal_point = [T.X_pix_(j+1)+1 T.Y_pix_(j+1)+1];
            %plot([start_point(1), goal_point(1)]',[start_point(2), goal_point(2)]','Marker','o')
            n = 1;
            t = linspace(0,1,n+2);  
            t = t(2:(end-1));       
            v = goal_point - start_point;
            x = start_point(1) + t*v(1);
            y = start_point(2) + t*v(2);
            v = w*v / norm(v);
            for i=1:n
                %plot([x(i)+v(2), x(i)-v(2)]',[y(i)-v(1), y(i)+v(1)]');
                c = improfile(Im_green,[x(i)+v(2), x(i)-v(2)]',[y(i)-v(1), y(i)+v(1)]');
            end
            distance = distance + pdist([start_point(1), start_point(2); goal_point(1), goal_point(2)],'euclidean');
            try
                f = fit([1:length(c)]',c-min(c),'gauss1');
            catch 
                f.c1 = 0;
            end
            if f.c1 ~= 0 && f.b1 > 0 && f.c1 < 2*w && f.b1 < 2*w
                %f.b1
                count = count+1;
                diameter_array(ID_count, count) = f.c1;
                dist_array(ID_count, count) = distance;
            end
        end
    end
end

edges =0.00000000001:0.05:1.000000000001; 
dia_hist = zeros(length(edges), 1);
for ii = 1:ID_count
    dist_array(ii, :) = dist_array(ii, :)/max(dist_array(ii, :));
    index = discretize(dist_array(ii, :), edges);
    diameter = diameter_array(ii, :);
    for jj = unique(index(~isnan(index)))
        zero_index = find(dia_hist(jj, :)==0);
        dia_hist(jj, zero_index : zero_index+length(diameter(jj == index))) = [diameter(jj == index), 0];
    end
end
dia_hist(dia_hist == 0) = NaN;
figure(12345);plot(mean(dia_hist, 2, 'omitnan')*2.35);hold on;
