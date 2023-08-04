function ClincFreq_Norm2MNI_T1_clinc(indat,inroi,outdir,outname,voxelsize,pth)
try
[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
% load([pth,filesep,'Norm2T1_segment.mat']);
load([pth,filesep,'Norm2T1_clincal.mat']);

pathspm = which('spm.m');
[patspm,namspm,extspm] = fileparts(pathspm);

MATLABBAT_seg = matlabbatch;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.anat{1} = indat;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.les{1} = inroi{1,1};
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.vox = [1,1,1]*voxelsize;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.bb = [-90 -126 -72;90 90 108];
spm_jobman('run_nogui',MATLABBAT_seg)

if ~isempty([patdat,filesep,'w',namdat,extdat])
    copyfile([patdat,filesep,'w',namdat,extdat],[outdir,filesep,'image_',outname,'_w',num2str(voxelsize),'_',namdat,extdat])
else
    
end
[patroi,namroi,extroi] = fileparts(inroi{1,1});
[v d] = Dynamic_read_dir_NIFTI([patroi,filesep,'ws',namroi,extdat]);
DAT = zeros(v.dim);
DAT(d>0.8) = 1;
DynamicBC_write_NIFTI(DAT,v,[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);
copyfile([patroi,filesep,'ws',namroi,extdat],[outdir,filesep,'smoothedROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);
catch
    
end
end