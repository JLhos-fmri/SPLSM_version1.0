function ClincFreq_WeiCentStat(outdirCent,OutWeiCentdir,SetUpPara,Para,AtlasTempU,dbmask)
IndirCent = [outdirCent,'sVolWeiCent',filesep];
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(IndirCent);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
dbmaskindn = find(dbmask==0);
% D = D>0;
for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    WeiCentDir{i} = [OutWeiCentdir,filesep,LabName,filesep];
    
    if (Para.WeiCent.V&&length(dir([WeiCentDir{i},'VoxelWise',filesep,'*.nii']))<2)||...
            (Para.WeiCent.R&&length(dir([WeiCentDir{i},'ROIWise',filesep,AtlasTempU{1,3},filesep,'Group*.csv']))<2)
        continue
    end
%     if Para.WeiCent.R&&length(dir([WeiCentDir{i},'ROIWise',filesep,AtlasTempU{1,3},filesep,'Group*.csv']))<2
%         continue
%     end
    
    %     mkdir(WeiCentDir{i});
    if Para.WeiCent.V
        VoxelDir = [WeiCentDir{i},'VoxelWise',filesep];
        %         mkdir(VoxelDir);
        if Para.WeiCent.StPermutation
            
            VoxelDirPerm = [WeiCentDir{i},'VoxelWise_Permutation',filesep];            
            VoxelDirPermcluster = [WeiCentDir{i},'VoxelWise_Permutationcluster',filesep];
            mkdir(VoxelDirPerm);
            mkdir(VoxelDirPermcluster);
            C = nchoosek(1:length(SetUpPara.ParaOut(i).datalab),2);
            for j = 1:size(C,1)
                IDtemp1 = SetUpPara.ParaOut(i).datalab{C(j,1)};
                IDtemp2 = SetUpPara.ParaOut(i).datalab{C(j,2)};
                DATtemp1 = full(D(:,IDtemp1));
                DATtemp2 = full(D(:,IDtemp2));
                DATtemp1(dbmaskindn,:) = 0;
                DATtemp2(dbmaskindn,:) = 0;
                NUM1 = length(IDtemp1);
                NUM2 = length(IDtemp2);
%                 [pval ClusNum] = Clinc_permtest(DATtemp1,DATtemp2,NUM1,NUM2,1,V(1));
                [pval ClusNum] = Clinc_permtest(DATtemp1,DATtemp2,NUM1,NUM2,1,V(1),VoxelDirPermcluster);                
%                 [pval ClusNum] = Clinc_permtest(DATtemp1,DATtemp2,NUM1,NUM2,1,V(1),VoxelDirPermcluster);
%                                 Clinc_permtest(DATtemp1,DATtemp2,NUM1,NUM2,Lab,V,VoxelDirPermcluster)
                Pvals = reshape(pval,V(1).dim);
                DynamicBC_write_NIFTI(Pvals,V(1),[VoxelDirPerm,['Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2))],'_PermP.nii'])
                save([VoxelDirPerm,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_CorrectClusNum.mat'],'ClusNum');
                delete([VoxelDirPermcluster,'*.mat'])
            end
        end
    end
    
    if Para.WeiCent.R
        ROIDir = [WeiCentDir{i},'ROIWise',filesep];
        
        if Para.WeiCent.StPermutation
            ROIDirPerm = [WeiCentDir{i},'ROIWise_Permutation',filesep];
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
                    DATTemp1 = G1dat.IND_involve1;
                    DATTemp2 = G2dat.IND_involve1;
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