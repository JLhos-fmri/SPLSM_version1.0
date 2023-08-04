function ClincFreq_CalClinWeiLocation(INDIR,OutWeiLocdir,SetUpPara,Para,VolLab,ClincWei)
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

for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    OWLDir{i} = [OutWeiLocdir,filesep,LabName,filesep,'VoxelWise',filesep];
    mkdir(OWLDir{i});
    
    for j = 1:length(SetUpPara.ParaOut(i).datalab)
        IDtemp = SetUpPara.ParaOut(i).datalab{j};
        GroupNum(j,1) = length(IDtemp);
        DATtemp = D(:,IDtemp);
        DATindex = find(sum(DATtemp,2)>0);
        for k = 1:GroupNum(j,1)
            vol = nnz(DATtemp(:,k))*Vsize^3/1000;
            indtemp = find(DATtemp(:,k));
            if VolLab
                VolWeiDATtemp(indtemp,k) = vol;
            else
                VolWeiDATtemp(indtemp,k) = ClincWei.Var(IDtemp(k));
            end
        end
        VolWeiDATindex_exist = VolWeiDATtemp(DATindex,:);
        parfor k = 1:length(DATindex)
            DATTEMP = VolWeiDATindex_exist(k,:);
            [mediaV(k),meanV(k),stdV(k)] = GiveClinc2Vol(DATTEMP);            
        end
        VolWeiOut_median = zeros(V(1).dim);
        VolWeiOut_median(DATindex) = mediaV;
        VolWeiOut_mean = zeros(V(1).dim);
        VolWeiOut_mean(DATindex) = meanV;
        VolWeiOut_std = zeros(V(1).dim);
        VolWeiOut_std(DATindex) = stdV;
        if VolLab
            DynamicBC_write_NIFTI(VolWeiOut_median,V(1),[OWLDir{i},'VolumeWeighted_',LabName,'_Group',sprintf('%03d',j),'_median.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_mean,V(1),[OWLDir{i},'VolumeWeighted_',LabName,'_Group',sprintf('%03d',j),'_mean.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_std,V(1),[OWLDir{i},'VolumeWeighted_',LabName,'_Group',sprintf('%03d',j),'_std.nii'])
            %
            DynamicBC_write_NIFTI(VolWeiOut_median.*Dmasked,V(1),[OWLDir{i},'IG_VolumeWeighted_',LabName,'_Group',sprintf('%03d',j),'_median.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_mean.*Dmasked,V(1),[OWLDir{i},'IG_VolumeWeighted_',LabName,'_Group',sprintf('%03d',j),'_mean.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_std.*Dmasked,V(1),[OWLDir{i},'IG_VolumeWeighted_',LabName,'_Group',sprintf('%03d',j),'_std.nii'])
        else            
            DynamicBC_write_NIFTI(VolWeiOut_median,V(1),[OWLDir{i},ClincWei.Name,'Weighted_',LabName,'_Group',sprintf('%03d',j),'_median.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_mean,V(1),[OWLDir{i},ClincWei.Name,'Weighted_',LabName,'_Group',sprintf('%03d',j),'_mean.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_std,V(1),[OWLDir{i},ClincWei.Name,'Weighted_',LabName,'_Group',sprintf('%03d',j),'_std.nii'])
            %
            DynamicBC_write_NIFTI(VolWeiOut_median.*Dmasked,V(1),[OWLDir{i},'IG_',ClincWei.Name,'Weighted_',LabName,'_Group',sprintf('%03d',j),'_median.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_mean.*Dmasked,V(1),[OWLDir{i},'IG_',ClincWei.Name,'Weighted_',LabName,'_Group',sprintf('%03d',j),'_mean.nii'])
            DynamicBC_write_NIFTI(VolWeiOut_std.*Dmasked,V(1),[OWLDir{i},'IG_',ClincWei.Name,'Weighted_',LabName,'_Group',sprintf('%03d',j),'_std.nii'])
        end
        clear mediaV meanV stdV
    end
end

end

function [mediaV,meanV,stdV] = GiveClinc2Vol(DATTEMP)

DATU = DATTEMP(DATTEMP>0);
mediaV = median(DATU);
meanV = mean(DATU);
stdV = std(DATU);

end