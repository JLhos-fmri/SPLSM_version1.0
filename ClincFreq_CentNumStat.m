function ClincFreq_CentNumStat(outdirCent,OutCentdir,SetUpPara,Para,AtlasTempU)

for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    CentDir{i} = [OutCentdir,filesep,LabName,filesep];
    if Para.CentNum.R
        if length(dir([CentDir{i},'ROIWise',filesep,AtlasTempU{1,3},filesep,'Group*.csv']))<2
            continue
        end
    end
    
    if Para.CentNum.R
        ROIDir = [CentDir{i},'ROIWise',filesep];
        
        if Para.CentNum.StChisquare            
            ROIDirChi2test = [CentDir{i},'ROIWise_Chi2test',filesep];
            mkdir(ROIDirChi2test);
            for iROI = 1:size(AtlasTempU,1)
                dirROI = [ROIDir,AtlasTempU{iROI,3},filesep];
                %
                ROIDirChi2testTemp = [ROIDirChi2test,AtlasTempU{iROI,3},filesep];
                mkdir(ROIDirChi2testTemp);
                %
                markedAtlas = dir([dirROI,'Group*.csv']);
                clear MAtlas MBoth MROI Nums
                for igroup = 1:length(markedAtlas)
                    aa = markedAtlas(igroup).name;
                    indsep = find(aa=='-');
                    Nums(igroup,1) = str2num(aa(indsep(end)+6:end-4));
                    MAtlas(igroup,:) = load([dirROI,markedAtlas(igroup).name]);
                end
                C = nchoosek(1:length(markedAtlas),2);
                
                MAtlas2 = repmat(Nums,1,size(MAtlas,2))-MAtlas;
                clear p_ANOVA_MAtlas Q_ANOVA_MAtlas p_two_MAtlas q_two_MAtlas or_two_MAtlas or_twolow_MAtlas or_twohigh_MAtlas
                for idat = 1:size(MAtlas,2)
                    x = [MAtlas(:,idat),MAtlas2(:,idat)];
                    [p_ANOVA_MAtlas(idat),Q_ANOVA_MAtlas(idat),p_two_MAtlas(idat,:),...
                        q_two_MAtlas(idat,:),or_two_MAtlas(idat,:),or_twolow_MAtlas(idat,:),...
                        or_twohigh_MAtlas(idat,:)] = ClincFreq_Chi2test(x,C);
                end
                [VROI DROI] = Dynamic_read_dir_NIFTI(AtlasTempU{iROI,1});
                
                if size(C,1)>1
                    csvwrite([ROIDirChi2testTemp,'p_ANOVA.csv'],p_ANOVA_MAtlas);
                    csvwrite([ROIDirChi2testTemp,'Q_ANOVA.csv'],Q_ANOVA_MAtlas);
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_ANOVA_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_ANOVA.nii']);
                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Q_ANOVA_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_ANOVA.nii']);                    
                end
                for iroi = 1:size(C,1)
                    %
                    csvwrite([ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.csv'],p_two_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_two_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.nii']);
                    csvwrite([ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.csv'],q_two_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = q_two_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.nii']);
                    
                    
                    csvwrite([ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.csv'],or_two_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_two_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.nii']);
                    
                    csvwrite([ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.csv'],or_twolow_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twolow_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.nii']);
                    
                    csvwrite([ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.csv'],or_twohigh_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twohigh_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'.nii']);
                end
            end
        end
        
        if Para.Heat.StPermutation
            ROIDirPerm = [CentDir{i},'ROIWise_Permutation',filesep];
            mkdir(ROIDirPerm);
            for iROI = 1:size(AtlasTempU,1)
                ROIdirss = [ROIDir,AtlasTempU{iROI,3},filesep];
                mats = dir([ROIdirss,'Single*.mat']);
                ROIDirPermU = [ROIDirPerm,AtlasTempU{iROI,3},filesep];
                mkdir(ROIDirPermU);
                
                [VROI DROI] = Dynamic_read_dir_NIFTI(AtlasTempU{iROI,1});
                
                C = nchoosek(1:length(mats),2);
                for iC = 1:size(C,1)
                    G1dat = load([ROIdirss,mats(C(iC,1)).name]);
                    G2dat = load([ROIdirss,mats(C(iC,2)).name]);
                    DATTemp1 = G1dat.INDinvolv1;
                    DATTemp2 = G2dat.INDinvolv1;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MAtlas = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);
                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedAtlas.csv'],Permpval_MAtlas);
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedAtlas.nii']);
                    
                    clear Permpval_MAtlas Permpval_MBoth Permpval_MROI
                end
            end
        end
    end
end

end