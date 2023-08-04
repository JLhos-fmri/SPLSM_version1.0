function ClincFreq_CentPointForVolume_pre(INDIR,outdirCent,Para,spmjobdir)
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;
inum = 1;
for isub = 1:length(namelist)
    Du = reshape(full(D(:,isub)),V(1).dim);
    Coninf = bwconncomp(Du);
    
    nums = Coninf.NumObjects;
    if nums>0
        for i = 1:nums
            Dattemp = zeros(size(Du));
            Dattemp(Coninf.PixelIdxList{i}) = 1;
            stats = regionprops3(Dattemp);
            Cent(i,:) = stats.Centroid;
            Vol(i,:) = stats.Volume;
        end
        
        Data_pos = zeros(size(Du));
        Data_vol = zeros(size(Du));
        for i = 1:nums
            Data_pos(round(Cent(i,2)),round(Cent(i,1)),round(Cent(i,3))) = 1;
            Data_vol(round(Cent(i,2)),round(Cent(i,1)),round(Cent(i,3))) = Vol(i);
        end
        
        CentOut = [round(Cent(:,2)),round(Cent(:,1)),round(Cent(:,3))];
        CentOutMNI = cor2mni(CentOut,V(1).mat);
        Outmat.CentOut = CentOut;
        Outmat.CentOutMNI = CentOutMNI;
        Outmat.Volume = Vol*abs(V(1,1).mat(1)*V(1).mat(2,2)*V(1).mat(3,3))/1000;% with mL as unit  cm3
        prename = namelist{isub};
        DynamicBC_write_NIFTI(Data_pos,V(1),[outdirCent,prename(1:end-4),'_Position.nii']);
        if Para.WeiCent.V||Para.WeiCent.R
            DynamicBC_write_NIFTI(Data_vol,V(1),[outdirCent,prename(1:end-4),'_Volume.nii']);
            
            namelistU{inum,1} = [outdirCent,prename(1:end-4),'_Volume.nii'];
            inum = inum+1;
        end
        save([outdirCent,prename(1:end-4),'_outinfor.mat'],'Outmat');
    else
        Outmat.CentOut = 0;
        Outmat.CentOutMNI = 0;
        Outmat.Volume = 0;
        
        prename = namelist{isub};
        save([outdirCent,prename(1:end-4),'_outinfor.mat'],'Outmat');
        Data_pos = zeros(V(1).dim);
        Data_vol = zeros(V(1).dim);
        
        DynamicBC_write_NIFTI(Data_pos,V(1),[outdirCent,prename(1:end-4),'_Position.nii']);
        if Para.WeiCent.V||Para.WeiCent.R
            DynamicBC_write_NIFTI(Data_vol,V(1),[outdirCent,prename(1:end-4),'_Volume.nii']);
            DynamicBC_write_NIFTI(Data_vol,V(1),[outdirCent,'s',prename(1:end-4),'_Volume.nii']);
        end
    end
end

if Para.WeiCent.V||Para.WeiCent.R
    load(fullfile(spmjobdir,'SmoothForVolume.mat'));
    spm_jobman('initcfg');   
%     matlabbatch
%     spm('init')
    for isub = 1:length(namelistU)        
        MATLABBAT = matlabbatch;
        MATLABBAT{1}.spm.spatial.smooth.data{1} = namelistU{isub,1};
        MATBAT{isub} = MATLABBAT;
    end
    parfor isub = 1:length(namelistU)
        MATBATU = MATBAT{isub};
        spm_jobman('run',MATBATU);
    end
end

OutdirVol1 = [outdirCent,'CentNum',filesep];
if Para.WeiCent.V||Para.WeiCent.R
    OutdirVol2 = [outdirCent,'VolWeiCent',filesep];
    OutdirVol3 = [outdirCent,'sVolWeiCent',filesep];
    mkdir(OutdirVol2);
    mkdir(OutdirVol3);
end
OutdirVol4 = [outdirCent,'Matfile',filesep];
mkdir(OutdirVol1);
mkdir(OutdirVol4);
movefile([outdirCent,'*_Position.nii'],OutdirVol1);
if Para.WeiCent.V||Para.WeiCent.R
    for i = 1:length(namelist)
        prename = namelist{i};
        file1 = [outdirCent,prename(1:end-4),'_Volume.nii'];
        file2 = [outdirCent,'s',prename(1:end-4),'_Volume.nii'];
        movefile(file1,OutdirVol2);
        movefile(file2,OutdirVol3);
    end
end
movefile([outdirCent,'*.mat'],OutdirVol4);
end