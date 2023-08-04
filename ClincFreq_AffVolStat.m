function ClincFreq_AffVolStat(INDIR,OutAffVoldir,SetUpPara,Para,AtlasTempU)

for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    AffVoldir{i} = [OutAffVoldir,filesep,LabName,filesep];
    
    
    if length(dir([AffVoldir{i},'ROIWise',filesep,AtlasTempU{1,3},filesep,'Group*AffVol.csv']))<2
        continue
    end
    
    if Para.AffVolume.R
        ROIDir = [AffVoldir{i},'ROIWise',filesep];
        
        if Para.AffVolume.StPermutation
            ROIDirPerm = [AffVoldir{i},'ROIWise_Permutation',filesep];
            mkdir(ROIDirPerm);
            for iROI = 1:size(AtlasTempU,1)
                ROIdirss = [ROIDir,AtlasTempU{iROI,3},filesep];
                mats = dir([ROIdirss,'Single*.mat']);
                ROIDirPermU = [ROIDirPerm,AtlasTempU{iROI,3},filesep];
                mkdir(ROIDirPermU);
                
                [VROI DROI] = Dynamic_read_dir_NIFTI(AtlasTempU{iROI,1});
                
                C = nchoosek(1:length(mats),2);
                for iC = 1:size(C,1)                    
                    indtemp = AtlasTempU{iROI,6};
                    
                    G1dat = load([ROIdirss,mats(C(iC,1)).name]);
                    G2dat = load([ROIdirss,mats(C(iC,2)).name]);
                    DATTemp1 = G1dat.INDinvolv1;
                    DATTemp2 = G2dat.INDinvolv1;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MAtlas = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedAtlas.csv'],Permpval_MAtlas);
                    DAT = zeros(VROI.dim);
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedAtlas.nii']);
                    
                    
                    DATTemp1 = G1dat.INDinvolv2;
                    DATTemp2 = G2dat.INDinvolv2;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MBoth = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedROI.csv'],Permpval_MBoth);
                    DAT = zeros(VROI.dim);
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MBoth(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedROI.nii']);
                    
                    DATTemp1 = G1dat.INDinvolv3;
                    DATTemp2 = G2dat.INDinvolv3;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MBoth = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedBoth.csv'],Permpval_MBoth);
                    DAT = zeros(VROI.dim);
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MBoth(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedBoth.nii']);
                    
                    %
                    DATTemp1 = G1dat.INDinvolv4;
                    DATTemp2 = G2dat.INDinvolv4;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_AffVol = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'AffVol.csv'],Permpval_AffVol);
                    DAT = zeros(VROI.dim);
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_AffVol(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_AffVol.nii']);
                    clear Permpval_MAtlas Permpval_MBoth Permpval_MROI Permpval_AffVol
                end
            end
        end
    end
end

end