clear all
close all

matlabStormPath = 'P:\MATLAB_Analysis\matlab-storm\';
display('Adding analysis code path...');
addpath('P:\MATLAB_Analysis\RNA_m6A_analysis\');

mainDirPath = 'L:\20180207 two color CB1 NCAM B2\DIV18 rabNCAM1-Cy3-647 B2-680\beads\split\';
InputFileName_l = [mainDirPath 'Conv_NCAM-Cy3-647_B2-680_0001-L_mlist.bin'];
InputFileName_r = [mainDirPath 'Conv_NCAM-Cy3-647_B2-680_0001-R_mlist.bin'];
OutputFileName = [mainDirPath 'tform.mat'];
OutputFileNamerev = [mainDirPath 'tformrev.mat'];

lmol = readbinfile_cmos(InputFileName_l);
figure(1)
plot(lmol.x,lmol.y,'ro');
title('original image');
xlim([0 256]);
ylim([0 256]);
axis equal;
hold on;

rmol = readbinfile_cmos(InputFileName_r);

plot(rmol.x,rmol.y,'b.');
xlim([0 256]);
ylim([0 256]);
axis equal;
hold off;

WhetherROI = questdlg('Do you want to select ROIs ?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
    selectROI;    
    lmol_match =[];
    rmol_match =[];
    for i = 1:length(roiList)
        x1 = roiList(i).rect(1);
        x2 = roiList(i).rect(1)+roiList(i).rect(3);
        y1 = roiList(i).rect(2);
        y2 = roiList(i).rect(2)+roiList(i).rect(4);
        lmol_select =[];
        rmol_select =[];
        for j = 1:length(lmol.x)
            if lmol.x(j)>x1 && lmol.x(j)<x2 && lmol.y(j)>y1 && lmol.y(j)<y2
                lmol_select = cat(1, lmol_select, [lmol.x(j) lmol.y(j)]);
            end
        end
        for j = 1:length(rmol.x)
            if rmol.x(j)>x1 && rmol.x(j)<x2 && rmol.y(j)>y1 && rmol.y(j)<y2
                rmol_select = cat(1, rmol_select, [rmol.x(j) rmol.y(j)]);
            end
        end
        lmol_match = cat(1, lmol_match, mean(lmol_select,1));
        rmol_match = cat(1, rmol_match, mean(rmol_select,1));
    end
    tform = cp2tform(lmol_match,rmol_match,'projective'); 
    tformrev = cp2tform(rmol_match,lmol_match,'projective'); 
    save(OutputFileName,'tform');
    save(OutputFileNamerev,'tformrev');
end


