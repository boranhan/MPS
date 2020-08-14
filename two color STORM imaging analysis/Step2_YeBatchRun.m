%Folder of original dax file, folder of bin file, bin file name, subfolder of beads tform in original dax file. 
   
for i = [1:6]
    FuncMasterRunFile('N:\STORM1_Data\20180122_TIAR_RBP\B1_TIAR-680_TIA1-647\', 'Splitted\', ['storm  add-647 b2-680_000' num2str(i)], 'Splitted\');
end

for i = [1:6]
    FuncMasterRunFile('N:\STORM1_Data\20180122_TIAR_RBP\B2_TIAR-680_eIF4A1-647\', 'Splitted\', ['STORM_647_000' num2str(i)], 'Splitted\');
end


for i = [1:5]
    FuncMasterRunFile('N:\STORM1_Data\20180122_TIAR_RBP\C1_G3BP1-680_YTHDF1-647\', 'Splitted\', ['STORM_647_000' num2str(i)], 'Splitted\');
end

for i = [2:9]
    FuncMasterRunFile('N:\STORM1_Data\20180125_NaAsO2release\A1_30+60_G3BP1-680_YTHDF1-647\', 'Splitted\', ['STORM_647_000' num2str(i)], 'Splitted\');
end
for i = [1:4]
    FuncMasterRunFile('N:\STORM1_Data\20180125_NaAsO2release\A3_30+60_G3BP1-680_YTHDF3-647\', 'Splitted\', ['STORM_647_000' num2str(i)], 'Splitted\');
end

for i = [1:7]
    FuncMasterRunFile('N:\STORM1_Data\20180125_NaAsO2release\A4_30+60_G3BP1-680_m6A-pAb-2-647\', 'Splitted\', ['STORM_647_000' num2str(i)], 'Splitted\');
end