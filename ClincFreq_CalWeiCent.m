function ClincFreq_CalWeiCent(outdirCent,OutWeiCentdir,SetUpPara,Para,AtlasTempU,dbmask)
% [V0,D0,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
% D0(isnan(D0)) = 0;
% D0(isinf(D0)) = 0;
% D0 = D0>0;
% 
% Dvol = sum(D0)';
% Dvole = sum(D0,2);
% IgV = Para.IgP;
% indexs = find(Dvole>=(size(D0,2)*IgV/100));
% Dmasked = zeros(V0.dim);
% Dmasked(indexs) = 1;
%%

IndirCent = [outdirCent,'sVolWeiCent',filesep];

[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(IndirCent);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
% D = D*abs(V(1).mat(1,1)*V(1).mat(2,2)*V(1).mat(3,3))/1000;


for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    WeiCentDir{i} = [OutWeiCentdir,filesep,LabName,filesep];
    mkdir(WeiCentDir{i});
    if Para.WeiCent.V
        VoxelDir = [WeiCentDir{i},'VoxelWise',filesep];
        mkdir(VoxelDir);
    end
    if Para.WeiCent.R
        ROIDir = [WeiCentDir{i},'ROIWise',filesep];
        mkdir(ROIDir);
        for iROI = 1:size(AtlasTempU,1)
            mkdir([ROIDir,AtlasTempU{iROI,3},filesep]);
        end
    end
    for j = 1:length(SetUpPara.ParaOut(i).datalab)
        IDtemp = SetUpPara.ParaOut(i).datalab{j};
%         GroupNum(j,1) = length(IDtemp);
        DATtemp = D(:,IDtemp);
        if Para.WeiCent.V
            DATtempV = reshape(full(sum(DATtemp,2)),V(1).dim).*dbmask;
            DynamicBC_write_NIFTI(DATtempV,V(1),[VoxelDir,filesep,LabName,'_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.nii'])
        end
        if Para.WeiCent.R
            for iROI = 1:size(AtlasTempU,1)
%                 PercentVol = str2num(Para.WeiCent.Rper)/100;
                indtab = AtlasTempU{iROI,6};
                for k = 1:length(IDtemp)
                    DATsingle = full(DATtemp(:,k));
%                     indinfor = AtlasTempU{iROI,4};
                    for l = 1:length(indtab)
                        IND_involve1(k,l) = sum(DATsingle(indtab{l}));
                    end
                end
                
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.csv'],IND_involve1);
                
                NamelistSave = namelist(IDtemp);
%                 DATe = zeros(V(1).dim);
                save([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.mat'],'NamelistSave','V','IND_involve1');
                clear NamtlistSave
                
                IND_total1 = sum(IND_involve1,1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.csv'],IND_total1);
                
                
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_total1(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'.nii']);
                
                clear INDinvolv1 INDinvolv2 INDinvolv3 IND_involve1 IND_involve2 IND_involve3 IND_total1 IND_total2 IND_total3 indtab
            end            
        end        
    end    
end
end