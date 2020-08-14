%Folder of original dax file, folder of bin file, bin file name, subfolder of beads tform in original dax file. 
   
ipath='P:\20190301 my10 647 b2 680\biolegend rab myh10 c-term DIV23\';
filename = dir(fullfile(strcat([ipath '\'],'storm*.dax')));
for i = 1:length(filename)
    FuncMasterRunFile_insight(ipath, 'split\', filename(i).name(1:end-4), 'split\');
end

% for i = 10:13
%     FuncMasterRunFile('H:\20180319 two color CB1 Src D-Kv\rtDIV16 CB1-D-647 B2-A-680\', 'split\', ['STORM_CB1-D-647 B2-A-680_00' num2str(i)], 'split\');
% end
% 
% for i = 1:5
%     FuncMasterRunFile('H:\20180319 two color CB1 Src D-Kv\rtDIV16 FGFR-D-647 B2-A-680 10min\', 'split\', ['STORM_ FGFR-K-647 CB1-A-680_000' num2str(i)], 'split\');
% end
% 
% % for i = 10:15
% %     FuncMasterRunFile('L:\20180312 two color YEEI Ntrk2 FGFR\msDIV2-20 4ul DM-GFP-647 B2-680\', 'split\', ['STORM_DM-GFP-647_B2-680_00' num2str(i)], 'split\');
% % end
% 
% for i = 1:7
%     FuncMasterRunFile('H:\20180319 two color CB1 Src D-Kv\rtDIV16 Ntrk2-D-647 B2-A-680 10min\', 'split\', ['STORM_ Ntrk2-D-647 CB1-A-680_000' num2str(i)], 'split\');
% end
% % 
% for i = 10:11
%     FuncMasterRunFile('L:\20180312 two color YEEI Ntrk2 FGFR\msDIV2-20 5ul Ntrk2-GFP-647 CB1-680\', 'split\', ['STORM_Ntrk2-647_CB1-680_00' num2str(i)], 'split\');
% end
% for i = 13
%     FuncMasterRunFile('L:\20180312 two color YEEI Ntrk2 FGFR\msDIV2-20 5ul Ntrk2-GFP-647 CB1-680\', 'split\', ['STORM_Ntrk2-647_CB1-680_00' num2str(i)], 'split\');
% end
% for i = 15
%     FuncMasterRunFile('L:\20180312 two color YEEI Ntrk2 FGFR\msDIV2-20 5ul Ntrk2-GFP-647 CB1-680\', 'split\', ['STORM_Ntrk2-647_CB1-680_00' num2str(i)], 'split\');
% end
% for i = 17:20
%     FuncMasterRunFile('L:\20180312 two color YEEI Ntrk2 FGFR\msDIV2-20 5ul Ntrk2-GFP-647 CB1-680\', 'split\', ['STORM_Ntrk2-647_CB1-680_00' num2str(i)], 'split\');
% end