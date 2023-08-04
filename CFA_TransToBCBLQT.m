function CFA_TransToBCBLQT
PG = uigetdir(pwd,'Normalized ROI Directory');
dir2 = [PG,filesep,'NormROI',filesep];
dirout = [PG,filesep,'NormROI_BCBLQT',filesep];
mkdir(dirout);
disp('translation start!')
[pth nam ext] = fileparts(which('CFA_TransToBCBLQT.m'));
refimg = fullfile(pth,'MNI152_T1_1mm_182218182.nii');
files = dir([dir2,'*.nii']);
for i = 1:length(files)
    dynamicBC_Reslice([dir2,files(i).name],[dirout,files(i).name],[1 1 1],0,refimg);
end
disp('Translation finished!')
end