function ClincFreq_Norm2MNI_EPI(indat,inroi,outdir,outname,voxelsize,pth)
try
[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
[v1 d1] = Dynamic_read_dir_NIFTI(indat);
[v2 d2] = Dynamic_read_dir_NIFTI(inroi{1});
d1(isnan(d1)) = 0;
d1(isinf(d1)) = 0;
D1 = zeros(v1.dim);
D1(d1>0) = 1;
d2(isnan(d2)) = 0;
d2(isinf(d2)) = 0;
D2 = zeros(v2.dim);
D2(d2>0) = 1;
D3 = (D1-D2)>0;
maskused = [patdat,filesep,'mask_',namdat,'.nii'];
DynamicBC_write_NIFTI(D3,v1,maskused);

load([pth,filesep,'Norm2EPI_est.mat']);
pathspm = which('spm.m');
[patspm,namspm,extspm] = fileparts(pathspm);

MATLABBAT_est = matlabbatch;
MATLABBAT_est{1}.spm.tools.oldnorm.est.subj.source{1} = indat;
MATLABBAT_est{1}.spm.tools.oldnorm.est.subj.wtsrc = [];
MATLABBAT_est{1}.spm.tools.oldnorm.est.subj.wtsrc{1} = maskused;
MATLABBAT_est{1}.spm.tools.oldnorm.est.eoptions.template{1} = fullfile(patspm,'toolbox','OldNorm','EPI.nii');

spm_jobman('run',MATLABBAT_est);

load([pth,filesep,'Norm2EPIandCT_writeImage.mat'])
MATLABBAT_writeimage = matlabbatch;
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.subj.matname{1} = fullfile(patdat,[namdat,'_sn.mat']);
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.subj.resample{1} = indat;
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.roptions.vox = [1,1,1]*voxelsize;
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.roptions.prefix = ['w',num2str(voxelsize)];
spm_jobman('run',MATLABBAT_writeimage);

load([pth,filesep,'Norm2EPIandCT_writeROI.mat'])
MATLABBAT_writeroi = matlabbatch;
MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.subj.matname{1} = fullfile(patdat,[namdat,'_sn.mat']);
for i = 1:length(inroi)
    MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.subj.resample{i,1} = inroi{i,1};
end
MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.roptions.vox = [1,1,1]*voxelsize;
MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.roptions.prefix = ['w',num2str(voxelsize)];
spm_jobman('run',MATLABBAT_writeroi);

copyfile([patdat,filesep,'w',num2str(voxelsize),namdat,extdat],[outdir,filesep,'image_',outname,'_w',num2str(voxelsize),'_',namdat,extdat]);
% copyfile([patroi,filesep,'w',num2str(voxelsize),namroi,extdat],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);
for i = 1:length(inroi)
    [patroi,namroi,extroi] = fileparts(inroi{i,1});
    copyfile([patroi,filesep,'w',num2str(voxelsize),namroi,extdat],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);
end
catch
    
end
end