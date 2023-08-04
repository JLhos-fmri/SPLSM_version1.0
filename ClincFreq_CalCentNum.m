function ClincFreq_CalCentNum(outdirCent,OutCentdir,SetUpPara,Para,AtlasTempU)
IndirCent = [outdirCent,'CentNum',filesep];

[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(IndirCent);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
% D = D>0;

% Dvol = sum(D)';
% Dvole = sum(D,2);
% IgV = Para.IgP;
% indexs = find(Dvole>=(size(D,2)*IgV/100));
% Dmasked = zeros(V.dim);
% Dmasked(indexs) = 1;

for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    CentDir{i} = [OutCentdir,filesep,LabName,filesep];
    mkdir(CentDir{i});
    
    if Para.CentNum.R
        ROIDir = [CentDir{i},'ROIWise',filesep];
        mkdir(ROIDir);
        for iROI = 1:size(AtlasTempU,1)
            mkdir([ROIDir,AtlasTempU{iROI,3},filesep]);
        end
    end
    for j = 1:length(SetUpPara.ParaOut(i).datalab)
        IDtemp = SetUpPara.ParaOut(i).datalab{j};
        GroupNum(j,1) = length(IDtemp);
        DATtemp = D(:,IDtemp);

        if Para.CentNum.R
            for iROI = 1:size(AtlasTempU,1)
                
                indtab = AtlasTempU{iROI,6};
                for k = 1:length(IDtemp)
                    DATsingle = full(DATtemp(:,k));
%                     indinfor = AtlasTempU{iROI,4};
                    for l = 1:length(indtab)
                        INDinvolv1(k,l) = nnz(DATsingle(indtab{l}));
                    end
                end
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.csv'],INDinvolv1);
                NamelistSave = namelist(IDtemp);
%                 DATe = zeros(V(1).dim);
                save([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.mat'],'NamelistSave','V','INDinvolv1');
                clear NamtlistSave
                
                IND_total1 = sum(INDinvolv1,1);
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