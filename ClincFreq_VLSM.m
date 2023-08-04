function ClincFreq_VLSM(INDIR,OutVLSM,SetUpPara,Para,dbmask,ParaVLSM,COVVLSM)
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;
Dvol = sum(D)';
Dvole = sum(D,2);
IgV = Para.IgP;
indexs = find(Dvole>=(size(D,2)*IgV/100));
if ~isempty(indexs)
    for i = 1:size(ParaVLSM,1)
        LabName = ParaVLSM{i,1};
        VALt = ParaVLSM{i,2};
        
        if Para.VLSM.COV
            for j = 1:size(COVVLSM,2)
                COVterm = [Dvol,COVVLSM(:,j)];
            end
        else
            COVterm = Dvol;
        end
        VLSMDir = [OutVLSM,filesep,LabName,filesep];
        mkdir(VLSMDir);
        VoxelDir = [VLSMDir,'VoxelWise',filesep];
        mkdir(VoxelDir);
        for j = 1:size(D,2)
            GROUP1{j,1} = 'G1';
            GROUP2{j,1} = 'G2';
        end
        Dind = D(indexs,:);
        t = zeros(length(indexs),1);
        p = ones(length(indexs),1)*0.5;
        Z = zeros(length(indexs),1);
        P = ones(length(indexs),1)*0.5;
        parfor j = 1:length(indexs)
            DAT = Dind(j,:);
            [t(j),p(j),Z(j),P(j)] = VLSMsubfunc(DAT,GROUP1,GROUP2,COVterm,VALt);
        end
        Tmap = zeros(V(1).dim);
        Pmap = zeros(V(1).dim);
        Zmap = zeros(V(1).dim);
        Tmap(indexs) = t;
        Pmap(indexs) = p;
        Zmap(indexs) = Z;
        DynamicBC_write_NIFTI(Tmap,V(1),[VoxelDir,filesep,'Tmap.nii']);
        DynamicBC_write_NIFTI(Pmap,V(1),[VoxelDir,filesep,'Pmap.nii']);
        DynamicBC_write_NIFTI(Zmap,V(1),[VoxelDir,filesep,'Zmap.nii']);
        save([VoxelDir,filesep,'ForPerm.mat'],'Dind','indexs','GROUP1','GROUP2','COVterm','VALt','V')
    end
else
    disp('no voxel exist in VLSM analysis')
end
end

