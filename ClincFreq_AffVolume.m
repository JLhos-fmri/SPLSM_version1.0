function ClincFreq_AffVolume(INDIR,OutAffVoldir,SetUpPara,Para,AtlasTempU)
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;

% Dvol = sum(D)';
% Dvole = sum(D,2);
% IgV = Para.IgP;
% indexs = find(Dvole>=(size(D,2)*IgV/100));
% Dmasked = zeros(V.dim);
% Dmasked(indexs) = 1;

for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    AffVDir{i} = [OutAffVoldir,filesep,LabName,filesep];
    mkdir(AffVDir{i});
    if Para.AffVolume.R
        ROIDir = [AffVDir{i},'ROIWise',filesep];
        mkdir(ROIDir);
        for iROI = 1:size(AtlasTempU,1)
            mkdir([ROIDir,AtlasTempU{iROI,3},filesep]);
        end
    end
    for j = 1:length(SetUpPara.ParaOut(i).datalab)
        IDtemp = SetUpPara.ParaOut(i).datalab{j};
        GroupNum(j,1) = length(IDtemp);
        DATtemp = D(:,IDtemp);
        if Para.AffVolume.R
            for iROI = 1:size(AtlasTempU,1)
%                 PercentVol = str2num(Para.Heat.Rper)/100;
                indtab = AtlasTempU{iROI,6};
                for k = 1:length(IDtemp)
                    DATsingle = full(DATtemp(:,k));
%                     indinfor = AtlasTempU{iROI,4};
                    for l = 1:length(indtab)
                        INDinvolv1(k,l) = nnz(DATsingle(indtab{l}))/nnz(DATsingle);
                        INDinvolv2(k,l) = nnz(DATsingle(indtab{l}))/length(indtab{l});
                        INDinvolv3(k,l) = max(INDinvolv1(k,l),INDinvolv2(k,l));
                        INDinvolv4(k,l) = nnz(DATsingle(indtab{l}))*abs(V(1).mat(1,1)*V(1).mat(2,2)*V(1).mat(3,3))/1000;
                    end
                end
%                 IND_involve1 = INDinvolv1>PercentVol;
%                 IND_involve2 = INDinvolv2>PercentVol;
%                 IND_involve3 = INDinvolv3>PercentVol;
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedROI.csv'],INDinvolv1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedAtlas.csv'],INDinvolv2);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedBoth.csv'],INDinvolv3);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-AffVol.csv'],INDinvolv4);
                NamelistSave = namelist(IDtemp);
%                 DATe = zeros(V(1).dim);
                save([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.mat'],'NamelistSave','V','INDinvolv1','INDinvolv2','INDinvolv3','INDinvolv4');
                clear NamtlistSave
%                 ROIdirtempsingle = [ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-ROImarked',filesep];
%                 mkdir(ROIdirtempsingle);
%                 for isub = 1:length(IDtemp)
%                     DATu = zeros(V(1).dim);
%                     for l = 1:length(indtab)
%                         DATu(indtab{l}) = IND_involve1(isub,l);
%                     end
%                     DynamicBC_write_NIFTI(DATu,V(1),[ROIdirtempsingle,namelist{isub}]);
%                 end
                %
%                 ROIdirtempsingle = [ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-Atlasmarked',filesep];
%                 mkdir(ROIdirtempsingle);
%                 for isub = 1:length(IDtemp)
%                     DATu = zeros(V(1).dim);
%                     for l = 1:length(indtab)
%                         DATu(indtab{l}) = IND_involve2(isub,l);
%                     end
%                     DynamicBC_write_NIFTI(DATu,V(1),[ROIdirtempsingle,namelist{isub}]);
%                 end
                %
%                 ROIdirtempsingle = [ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-Bothmarked',filesep];
%                 mkdir(ROIdirtempsingle);
%                 for isub = 1:length(IDtemp)
%                     DATu = zeros(V(1).dim);
%                     for l = 1:length(indtab)
%                         DATu(indtab{l}) = IND_involve3(isub,l);
%                     end
%                     DynamicBC_write_NIFTI(DATu,V(1),[ROIdirtempsingle,namelist{isub}]);
%                 end
                %
                
                
                IND_totalV1 = sum(INDinvolv1,1);
                IND_totalV2 = sum(INDinvolv2,1);
                IND_totalV3 = sum(INDinvolv3,1);
                IND_totalV4 = sum(INDinvolv4,1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedROI.csv'],IND_totalV1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedAtlas.csv'],IND_totalV2);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedBoth.csv'],IND_totalV3);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-AffVol.csv'],IND_totalV4);
                
                
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_totalV1(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-ROImarked.nii']);
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_totalV2(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-Atlasmarked.nii']);
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_totalV3(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-Bothmarked.nii']);
                
                
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_totalV4(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-AffVol.nii']);
                
                clear INDinvolv1 INDinvolv2 INDinvolv3 INDinvolv4 IND_involve1 IND_involve2 IND_involve3 IND_involve4 IND_total1 IND_total2 IND_total3 IND_total4 indtab
            end
            
        end
        
    end
    
end

end