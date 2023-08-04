function ClincFreq_ClinWeiLocationStat(INDIR,OutWeiLocdir,SetUpPara,Para,VolLab,ClincWei)
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;


Dvol = sum(D)';
Dvole = sum(D,2);
IgV = Para.IgP;
indexs = find(Dvole>=(size(D,2)*IgV/100));
Dmasked = zeros(V.dim);
Dmasked(indexs) = 1;

Vsize = SetUpPara.OutVsize;

for i = 2:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    OWLDir{i} = [OutWeiLocdir,filesep,LabName,filesep,'VoxelWise_Permutation',filesep];
    OWLclus{i} = [OutWeiLocdir,filesep,LabName,filesep,'VoxelWise_PermutationCluster',filesep];
    mkdir(OWLDir{i});
    mkdir(OWLclus{i});
    
    datlabs = SetUpPara.ParaOut(i).datalab;
    C = nchoosek(1:length(datlabs),2);
    
    
    for j = 1:size(C,1)
        clear Pval_media Pval_mean Pval_std
        IDtemp1 = SetUpPara.ParaOut(i).datalab{C(j,1)};
        IDtemp2 = SetUpPara.ParaOut(i).datalab{C(j,2)};
        
        DATtemp1 = D(:,IDtemp1);
        DATtemp2 = D(:,IDtemp2);
        
        GroupNum1 = length(IDtemp1);
        GroupNum2 = length(IDtemp2);
        
        
        DATindex1 = find(sum(DATtemp1,2)>0);
        DATindex2 = find(sum(DATtemp2,2)>0);  
        DATIND = unique([DATindex1;DATindex2]);
        VolWeiDATtemp1Full = zeros(size(D,1),GroupNum1);
%         VolWeiDATtemp1Full(I)
        VolWeiDATtemp2Full = zeros(size(D,1),GroupNum2);
        for k = 1:GroupNum1
            vol = nnz(DATtemp1(:,k))*Vsize^3/1000;
            indtemp1 = find(DATtemp1(:,k));
            if VolLab
                VolWeiDATtemp1Full(indtemp1,k) = vol;
            else
                VolWeiDATtemp1Full(indtemp1,k) = ClincWei.Var(IDtemp1(k));
            end
        end
        for k = 1:GroupNum2
            vol = nnz(DATtemp2(:,k))*Vsize^3/1000;
            indtemp2 = find(DATtemp2(:,k));
            if VolLab
                VolWeiDATtemp2Full(indtemp2,k) = vol;
            else
                VolWeiDATtemp2Full(indtemp2,k) = ClincWei.Var(IDtemp2(k));
            end
        end
        clear mediaV meanV stdV
        VolWeiDATindex1_exist = VolWeiDATtemp1Full(DATIND,:);
        parfor k = 1:length(DATIND)
            DATTEMP = VolWeiDATindex1_exist(k,:);
            [mediaV1(k),meanV1(k),stdV1(k)] = GiveClinc2Vol(DATTEMP);            
        end        
        VolWeiDATindex2_exist = VolWeiDATtemp2Full(DATIND,:);
        parfor k = 1:length(DATIND)
            DATTEMP = VolWeiDATindex2_exist(k,:);
            [mediaV2(k),meanV2(k),stdV2(k)] = GiveClinc2Vol(DATTEMP);            
        end
        %% make permutation distribution
%%
%         VolWeiDATall = [VolWeiDATtemp1Full,VolWeiDATtemp2Full];
%         
%         clear Perm_mediaV1P Perm_mediaV2P Perm_meanV1P Perm_meanV2P Perm_stdV1P Perm_stdV2P
%         for iperm = 1:1000
%             PermInd = randperm(GroupNum1+GroupNum2);
%             PermInd1 = PermInd(1:GroupNum1);
%             PermInd2 = PermInd(GroupNum1+1:GroupNum1+GroupNum2);
%             
%             VolWeiDATindex1_exist = VolWeiDATall(DATIND,PermInd1);
%             parfor k = 1:length(DATIND)
%                 DATTEMP = VolWeiDATindex1_exist(k,:);
%                 [mediaV1P(k),meanV1P(k),stdV1P(k)] = GiveClinc2Vol(DATTEMP);
%             end
%             VolWeiDATindex2_exist = VolWeiDATall(DATIND,PermInd2);
%             parfor k = 1:length(DATIND)
%                 DATTEMP = VolWeiDATindex2_exist(k,:);
%                 [mediaV2P(k),meanV2P(k),stdV2P(k)] = GiveClinc2Vol(DATTEMP);
%             end
%             Perm_mediaV1P(iperm,:) = mediaV1P;
%             Perm_mediaV2P(iperm,:) = mediaV2P;
%             Perm_meanV1P(iperm,:) = meanV1P;
%             Perm_meanV2P(iperm,:) = meanV2P;
%             Perm_stdV1P(iperm,:) = stdV1P;
%             Perm_stdV2P(iperm,:) = stdV2P;
%         end
%         Perm_mediaV1P(isnan(Perm_mediaV1P)) = 0;
%         Perm_mediaV2P(isnan(Perm_mediaV2P)) = 0;
%         deltMediaV = Perm_mediaV1P-Perm_mediaV2P;
%         Perm_meanV1P(isnan(Perm_meanV1P)) = 0;
%         Perm_meanV2P(isnan(Perm_meanV2P)) = 0;
%         deltMeanV = Perm_meanV1P-Perm_meanV2P;
%         Perm_stdV1P(isnan(Perm_stdV1P)) = 0;
%         Perm_stdV2P(isnan(Perm_stdV2P)) = 0;
%         deltstdV = Perm_stdV1P-Perm_stdV2P;
%         
%         mediaV1(isnan(mediaV1)) = 0;
%         mediaV2(isnan(mediaV2)) = 0;
%         DeltMedia = mediaV1-mediaV2;
%         meanV1(isnan(meanV1)) = 0;
%         meanV2(isnan(meanV2)) = 0;
%         DeltMean = meanV1-meanV2;
%         stdV1(isnan(stdV1)) = 0;
%         stdV2(isnan(stdV2)) = 0;
%         DeltStd = stdV1-stdV2;
%         
%         
% %         origDtemp = origd(CUTpiece{icut});
%         [mu sig muci sigci] = normfit(deltMediaV);
%         Pval_media = normcdf(DeltMedia,mu,sig);
%         indexs = find(mu==0&sig==0);
%         Pval_media(indexs) = 0.5;
%         [mu sig muci sigci] = normfit(deltMeanV);
%         Pval_mean = normcdf(DeltMean,mu,sig);
%         indexs = find(mu==0&sig==0);
%         Pval_mean(indexs) = 0.5;
%         [mu sig muci sigci] = normfit(deltstdV);
%         Pval_std = normcdf(DeltStd,mu,sig);
%         indexs = find(mu==0&sig==0);
%         Pval_std(indexs) = 0.5;
%%        update 2022.11.22 use segment for small RAM.

        VolWeiDATall = [VolWeiDATtemp1Full,VolWeiDATtemp2Full];
        DATINDNUM = length(DATIND);
        segs = 20000;
        SegNum = floor(DATINDNUM/segs);
        clear segind segind3
        for iseg = 1:SegNum
            segind{iseg} = DATIND(1+segs*(iseg-1):segs*iseg);
            segind3{iseg} = 1+segs*(iseg-1):segs*iseg;
        end
        if SegNum*segs<DATINDNUM
            segind{SegNum+1} = DATIND(segs*SegNum+1:DATINDNUM);
            segind3{SegNum+1} = segs*SegNum+1:DATINDNUM;
        end
        mediaV1(isnan(mediaV1)) = 0;
        mediaV2(isnan(mediaV2)) = 0;
        DeltMedia = mediaV1-mediaV2;
        meanV1(isnan(meanV1)) = 0;
        meanV2(isnan(meanV2)) = 0;
        DeltMean = meanV1-meanV2;
        stdV1(isnan(stdV1)) = 0;
        stdV2(isnan(stdV2)) = 0;
        DeltStd = stdV1-stdV2;
        permnums = 1000;
        save([OWLclus{i},filesep,'permnums.mat'],'permnums');
        mkdir([OWLclus{i},filesep,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),filesep]);
        for iseg = 1:length(segind)
            clear Perm_mediaV1P Perm_mediaV2P Perm_meanV1P Perm_meanV2P Perm_stdV1P Perm_stdV2P
            for iperm = 1:permnums
                PermInd = randperm(GroupNum1+GroupNum2);
                PermInd1 = PermInd(1:GroupNum1);
                PermInd2 = PermInd(GroupNum1+1:GroupNum1+GroupNum2);
                
                VolWeiDATindex1_exist = VolWeiDATall(segind{iseg},PermInd1);
                parfor k = 1:length(segind{iseg})
                    DATTEMP = VolWeiDATindex1_exist(k,:);
                    [mediaV1P(k),meanV1P(k),stdV1P(k)] = GiveClinc2Vol(DATTEMP);
                end
                VolWeiDATindex2_exist = VolWeiDATall(segind{iseg},PermInd2);
                parfor k = 1:length(segind{iseg})
                    DATTEMP = VolWeiDATindex2_exist(k,:);
                    [mediaV2P(k),meanV2P(k),stdV2P(k)] = GiveClinc2Vol(DATTEMP);
                end
                Perm_mediaV1P(iperm,:) = mediaV1P;
                Perm_mediaV2P(iperm,:) = mediaV2P;
                Perm_meanV1P(iperm,:) = meanV1P;
                Perm_meanV2P(iperm,:) = meanV2P;
                Perm_stdV1P(iperm,:) = stdV1P;
                Perm_stdV2P(iperm,:) = stdV2P;
                clear mediaV1P meanV1P stdV1P mediaV2P meanV2P stdV2P
            end
            Perm_mediaV1P(isnan(Perm_mediaV1P)) = 0;
            Perm_mediaV2P(isnan(Perm_mediaV2P)) = 0;
            deltMediaV = Perm_mediaV1P-Perm_mediaV2P;
            Perm_meanV1P(isnan(Perm_meanV1P)) = 0;
            Perm_meanV2P(isnan(Perm_meanV2P)) = 0;
            deltMeanV = Perm_meanV1P-Perm_meanV2P;
            Perm_stdV1P(isnan(Perm_stdV1P)) = 0;
            Perm_stdV2P(isnan(Perm_stdV2P)) = 0;
            deltstdV = Perm_stdV1P-Perm_stdV2P;     
            
            %         origDtemp = origd(CUTpiece{icut});
            [mu sig muci sigci] = normfit(deltMediaV);
            Pval_media0 = normcdf(DeltMedia(segind3{iseg}),mu,sig);
            indexs = find(mu==0&sig==0);
            Pval_media0(indexs) = 0.5;
            [mu sig muci sigci] = normfit(deltMeanV);
            Pval_mean0 = normcdf(DeltMean(segind3{iseg}),mu,sig);
            indexs = find(mu==0&sig==0);
            Pval_mean0(indexs) = 0.5;
            [mu sig muci sigci] = normfit(deltstdV);
            Pval_std0 = normcdf(DeltStd(segind3{iseg}),mu,sig);
            indexs = find(mu==0&sig==0);
            Pval_std0(indexs) = 0.5;
            
            Pval_media(segind3{iseg}) = Pval_media0;
            Pval_mean(segind3{iseg}) = Pval_mean0;
            Pval_std(segind3{iseg}) = Pval_std0;
            
            for iperm = 1:1000
                [mu sig muci sigci] = normfit(deltMediaV);
                Pval_media0 = normcdf(deltMediaV(iperm,:),mu,sig);
                indexs = find(mu==0&sig==0);
                Pval_media0(indexs) = 0.5;
                [mu sig muci sigci] = normfit(deltMeanV);
                Pval_mean0 = normcdf(deltMediaV(iperm,:),mu,sig);
                indexs = find(mu==0&sig==0);
                Pval_mean0(indexs) = 0.5;
                [mu sig muci sigci] = normfit(deltstdV);
                Pval_std0 = normcdf(deltMediaV(iperm,:),mu,sig);
                indexs = find(mu==0&sig==0);
                Pval_std0(indexs) = 0.5;
                Pval_media_perm(iperm,:) = Pval_media0;
                Pval_mean_perm(iperm,:) = Pval_mean0;
                Pval_std_perm(iperm,:) = Pval_std0;
            end
            save([OWLclus{i},filesep,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),filesep,'Piece',sprintf('%05d',iseg),'.mat'],'Pval_media_perm','Pval_mean_perm','Pval_std_perm','segind3','iseg','DATIND','V','Dmasked');
            clear Pval_media_perm Pval_mean_perm Pval_std_perm
        end
%%
        VolWeiOut_median = zeros(V(1).dim);
        VolWeiOut_median(DATIND) = Pval_media;
        VolWeiOut_mean = zeros(V(1).dim);
        VolWeiOut_mean(DATIND) = Pval_mean;
        VolWeiOut_std = zeros(V(1).dim);
        VolWeiOut_std(DATIND) = Pval_std;
        DynamicBC_write_NIFTI(VolWeiOut_median,V(1),[OWLDir{i},'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_median_PermP.nii'])
        DynamicBC_write_NIFTI(VolWeiOut_mean,V(1),[OWLDir{i},'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_mean_PermP.nii'])
        DynamicBC_write_NIFTI(VolWeiOut_std,V(1),[OWLDir{i},'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_std_PermP.nii'])

        
        DynamicBC_write_NIFTI(VolWeiOut_median.*Dmasked,V(1),[OWLDir{i},'IG_Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_median_PermP.nii'])
        DynamicBC_write_NIFTI(VolWeiOut_mean.*Dmasked,V(1),[OWLDir{i},'IG_Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_mean_PermP.nii'])
        DynamicBC_write_NIFTI(VolWeiOut_std.*Dmasked,V(1),[OWLDir{i},'IG_Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_std_PermP.nii'])

        ExtClustInfo_CWL(OWLclus{i},C,j,Dmasked,OWLDir{i});
    end
end

end

function [mediaV,meanV,stdV] = GiveClinc2Vol(DATTEMP)

DATU = DATTEMP(DATTEMP>0);
mediaV = median(DATU);
meanV = mean(DATU);
stdV = std(DATU);

end