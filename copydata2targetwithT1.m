function copydata2targetwithT1(PGout,ROIdirs,IMGdirs,T1dirs,namelist)
for i = 1:length(namelist)
    [pth1 fil1 ext1] = fileparts(IMGdirs{i});
    [pth2 fil2 ext2] = fileparts(ROIdirs{i});
    [pth3 fil3 ext3] = fileparts(T1dirs{i});
    Outdirs = [PGout,filesep,namelist{i},filesep];
    mkdir(Outdirs);
    if strcmp(ext1,'.nii')
        copyfile(IMGdirs{i},[Outdirs,'OtherMode.nii']);
    elseif strcmp(ext1,'.hdr')
        [v d] = Dynamic_read_dir_NIFTI(IMGdirs{i});
        D = reshape(d,v.dim);
        DynamicBC_write_NIFTI(D,v,[Outdirs,'OtherMode.nii']);
    elseif strcmp(ext1,'.gz')
        gunzip(IMGdirs{i});
        temps = IMGdirs{i};
        copyfile(temps(1:end-3),[Outdirs,'OtherMode.nii']);
        delete(temps(1:end-3))
    end
    
    if strcmp(ext2,'.nii')        
        copyfile(ROIdirs{i},[Outdirs,'ROI.nii']);
    elseif strcmp(ext2,'.hdr')        
        [v d] = Dynamic_read_dir_NIFTI(ROIdirs{i});
        D = reshape(d,v.dim);
        DynamicBC_write_NIFTI(D,v,[Outdirs,'ROI.nii']);
    elseif strcmp(ext2,'.gz')
        gunzip(ROIdirs{i});
        temps = ROIdirs{i};
        copyfile(temps(1:end-3),[Outdirs,'ROI.nii']);        
        delete(temps(1:end-3))
    end
    
    
    if strcmp(ext3,'.nii')        
        copyfile(T1dirs{i},[Outdirs,'T1.nii']);
    elseif strcmp(ext3,'.hdr')        
        [v d] = Dynamic_read_dir_NIFTI(T1dirs{i});
        D = reshape(d,v.dim);
        DynamicBC_write_NIFTI(D,v,[Outdirs,'T1.nii']);
    elseif strcmp(ext3,'.gz')
        gunzip(T1dirs{i});
        temps = T1dirs{i};
        copyfile(temps(1:end-3),[Outdirs,'T1.nii']);        
        delete(temps(1:end-3))
    end
end
end