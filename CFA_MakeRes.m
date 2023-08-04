function CFA_MakeRes
if ~isempty(dir('SetUpPara.mat'))
    load('SetUpPara.mat');
    load('Setup_TempAndPara.mat');
else
    [fil PG] = uigetfile('*.mat',pwd,'SetUpPara.mat');
    load(fullfile(PG,fil));
    load(fullfile(PG,'Setup_TempAndPara.mat'));
    cd(PG);
end

outdir = [SetUpPara.Outdir,filesep];
if isempty(dir(outdir))
    outdir = [pwd,filesep];
end
Atlasdirs = [outdir,filesep,'ReslicedAtlas',filesep];
try
    [Vatl Datl] = Dynamic_read_dir_NIFTI(Atlasdirs);
catch
    
end
for i = 2:length(SetUpPara.ParaOut)
    if Para.Heat.V
        if Para.Heat.StChisquare
            Voxdir = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise',filesep];
            files = dir([Voxdir,'Freq*.nii']);
            if length(files)>=2
                Voxdir1 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test',filesep];
                Voxdir1out05 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P05',filesep];
                mkdir(Voxdir1out05);
                Voxdir1out005 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P005',filesep];
                mkdir(Voxdir1out005);
                Voxdir1out01 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P01',filesep];
                mkdir(Voxdir1out01);
                Voxdir1out001 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P001',filesep];
                mkdir(Voxdir1out001);
                for j = 1:length(files)
                    [V(j) D(:,j)] = Dynamic_read_dir_NIFTI([Voxdir,files(j).name]);
                end
                D(isnan(D)) = 0;
                D(isinf(D)) = 0;
                if length(files)>2
                    [VAq DAq] = Dynamic_read_dir_NIFTI([Voxdir1,'ANOVA_q.nii']);
                    [VAp DAp] = Dynamic_read_dir_NIFTI([Voxdir1,'ANOVA_p.nii']);
                    DAq(isnan(DAq)) = 0; DAq(isinf(DAq)) = 0;
                    DAq(isnan(DAp)) = 0; DAq(isinf(DAp)) = 0;
                    THR1 = DAp<0.05&DAp>0;
                    THR2 = DAp<0.01&DAp>0;
                    THR3 = DAp<0.005&DAp>0;
                    THR4 = DAp<0.001&DAp>0;
                    DynamicBC_write_NIFTI(reshape(THR1,VAq.dim),VAq,[Voxdir1out05,'ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR2,VAq.dim),VAq,[Voxdir1out01,'ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR3,VAq.dim),VAq,[Voxdir1out005,'ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR4,VAq.dim),VAq,[Voxdir1out001,'ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR1.*DAq,VAq.dim),VAq,[Voxdir1out05,'ANOVA_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(THR2.*DAq,VAq.dim),VAq,[Voxdir1out01,'ANOVA_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(THR3.*DAq,VAq.dim),VAq,[Voxdir1out005,'ANOVA_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(THR4.*DAq,VAq.dim),VAq,[Voxdir1out001,'ANOVA_Q_masked.nii']);
                end
                C = nchoosek(1:length(files),2);
                for j = 1:size(C,1)
                    DIF = D(:,C(j,1))-D(:,C(j,2));
                    LABTEMP = ['Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2))];
                    [VQ DQ] = Dynamic_read_dir_NIFTI([Voxdir1,LABTEMP,'_q.nii']);
                    [VP DP] = Dynamic_read_dir_NIFTI([Voxdir1,LABTEMP,'_p.nii']);
                    [VOR DOR] = Dynamic_read_dir_NIFTI([Voxdir1,LABTEMP,'_OR.nii']);
                    [VORcihigh DORcihigh] = Dynamic_read_dir_NIFTI([Voxdir1,LABTEMP,'_ORcihigh.nii']);
                    [VORcilow DORcilow] = Dynamic_read_dir_NIFTI([Voxdir1,LABTEMP,'_ORcilow.nii']);
                    THR1P = DP<0.05&DP>0&DIF>0;
                    THR1N = DP<0.05&DP>0&DIF<0;
%                     THR1 = THR1P-THR1N;
                    THR1 = THR1P+THR1N;
                    DynamicBC_write_NIFTI(reshape(THR1.*DP,VQ.dim),VQ,[Voxdir1out05,LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR1.*DQ,VQ.dim),VQ,[Voxdir1out05,LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR1).*DOR,VQ.dim),VQ,[Voxdir1out05,LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR1).*DORcihigh,VQ.dim),VQ,[Voxdir1out05,LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR1).*DORcilow,VQ.dim),VQ,[Voxdir1out05,LABTEMP,'_ORcilow_masked.nii']);
                    
                    
                    THR2P = DP<0.05&DP>0&DIF>0;
                    THR2N = DP<0.05&DP>0&DIF<0;
                    THR2 = THR2P-THR2N;
                    DynamicBC_write_NIFTI(reshape(THR2.*DP,VQ.dim),VQ,[Voxdir1out01,LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR2.*DQ,VQ.dim),VQ,[Voxdir1out01,LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR2).*DOR,VQ.dim),VQ,[Voxdir1out01,LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR2).*DORcihigh,VQ.dim),VQ,[Voxdir1out01,LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR2).*DORcilow,VQ.dim),VQ,[Voxdir1out01,LABTEMP,'_ORcilow_masked.nii']);
                    
                    
                    THR3P = DP<0.05&DP>0&DIF>0;
                    THR3N = DP<0.05&DP>0&DIF<0;
                    THR3 = THR3P-THR3N;
                    DynamicBC_write_NIFTI(reshape(THR3.*DP,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR3.*DQ,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR3).*DOR,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR3).*DORcihigh,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR3).*DORcilow,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_ORcilow_masked.nii']);
                    
                    
                    THR4P = DP<0.05&DP>0&DIF>0;
                    THR4N = DP<0.05&DP>0&DIF<0;
                    THR4 = THR4P-THR4N;
                    DynamicBC_write_NIFTI(reshape(THR4.*DP,VQ.dim),VQ,[Voxdir1out001,LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR4.*DQ,VQ.dim),VQ,[Voxdir1out001,LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR4).*DOR,VQ.dim),VQ,[Voxdir1out001,LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR4).*DORcihigh,VQ.dim),VQ,[Voxdir1out001,LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR4).*DORcilow,VQ.dim),VQ,[Voxdir1out001,LABTEMP,'_ORcilow_masked.nii']);
                end
            end
            %%
            if length(files)>=2
                Voxdir1 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test',filesep];
                Voxdir1out05 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P05',filesep];
%                 mkdir(Voxdir1out05);
                Voxdir1out005 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P005',filesep];
%                 mkdir(Voxdir1out005);
                Voxdir1out01 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P01',filesep];
%                 mkdir(Voxdir1out01);
                Voxdir1out001 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Chi2test_P001',filesep];
%                 mkdir(Voxdir1out001);
                for j = 1:length(files)
                    [V(j) D(:,j)] = Dynamic_read_dir_NIFTI([Voxdir,'IG_',files(j).name]);
                end
                D(isnan(D)) = 0;
                D(isinf(D)) = 0;
                if length(files)>2
                    [VAq DAq] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_ANOVA_q.nii']);
                    [VAp DAp] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_ANOVA_p.nii']);
                    DAq(isnan(DAq)) = 0; DAq(isinf(DAq)) = 0;
                    DAq(isnan(DAp)) = 0; DAq(isinf(DAp)) = 0;
                    THR1 = DAp<0.05&DAp>0;
                    THR2 = DAp<0.01&DAp>0;
                    THR3 = DAp<0.005&DAp>0;
                    THR4 = DAp<0.001&DAp>0;
                    DynamicBC_write_NIFTI(reshape(THR1,VAq.dim),VAq,[Voxdir1out05,'IG_ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR2,VAq.dim),VAq,[Voxdir1out01,'IG_ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR3,VAq.dim),VAq,[Voxdir1out005,'IG_ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR4,VAq.dim),VAq,[Voxdir1out001,'IG_ANOVA_P_mask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR1.*DAq,VAq.dim),VAq,[Voxdir1out05,'IG_ANOVA_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(THR2.*DAq,VAq.dim),VAq,[Voxdir1out01,'IG_ANOVA_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(THR3.*DAq,VAq.dim),VAq,[Voxdir1out005,'IG_ANOVA_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(THR4.*DAq,VAq.dim),VAq,[Voxdir1out001,'IG_ANOVA_Q_masked.nii']);
                end
                C = nchoosek(1:length(files),2);
                for j = 1:size(C,1)
                    DIF = D(:,C(j,1))-D(:,C(j,2));
                    LABTEMP = ['Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2))];
                    [VQ DQ] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_',LABTEMP,'_q.nii']);
                    [VP DP] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_',LABTEMP,'_p.nii']);
                    [VOR DOR] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_',LABTEMP,'_OR.nii']);
                    [VORcihigh DORcihigh] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_',LABTEMP,'_ORcihigh.nii']);
                    [VORcilow DORcilow] = Dynamic_read_dir_NIFTI([Voxdir1,'IG_',LABTEMP,'_ORcilow.nii']);
                    THR1P = DP<0.05&DP>0&DIF>0;
                    THR1N = DP<0.05&DP>0&DIF<0;
%                     THR1 = THR1P-THR1N;
                    THR1 = THR1P+THR1N;
                    DynamicBC_write_NIFTI(reshape(THR1.*DP,VQ.dim),VQ,[Voxdir1out05,'IG_',LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR1.*DQ,VQ.dim),VQ,[Voxdir1out05,'IG_',LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR1).*DOR,VQ.dim),VQ,[Voxdir1out05,'IG_',LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR1).*DORcihigh,VQ.dim),VQ,[Voxdir1out05,'IG_',LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR1).*DORcilow,VQ.dim),VQ,[Voxdir1out05,'IG_',LABTEMP,'_ORcilow_masked.nii']);
                    
                    
                    THR2P = DP<0.05&DP>0&DIF>0;
                    THR2N = DP<0.05&DP>0&DIF<0;
                    THR2 = THR2P-THR2N;
                    DynamicBC_write_NIFTI(reshape(THR2.*DP,VQ.dim),VQ,[Voxdir1out01,'IG_',LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR2.*DQ,VQ.dim),VQ,[Voxdir1out01,'IG_',LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR2).*DOR,VQ.dim),VQ,[Voxdir1out01,'IG_',LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR2).*DORcihigh,VQ.dim),VQ,[Voxdir1out01,'IG_',LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR2).*DORcilow,VQ.dim),VQ,[Voxdir1out01,'IG_',LABTEMP,'_ORcilow_masked.nii']);
                    
                    
                    THR3P = DP<0.05&DP>0&DIF>0;
                    THR3N = DP<0.05&DP>0&DIF<0;
                    THR3 = THR3P-THR3N;
                    DynamicBC_write_NIFTI(reshape(THR3.*DP,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR3.*DQ,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR3).*DOR,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR3).*DORcihigh,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR3).*DORcilow,VQ.dim),VQ,[Voxdir1out005,LABTEMP,'_ORcilow_masked.nii']);
                    
                    
                    THR4P = DP<0.05&DP>0&DIF>0;
                    THR4N = DP<0.05&DP>0&DIF<0;
                    THR4 = THR4P-THR4N;
                    DynamicBC_write_NIFTI(reshape(THR4.*DP,VQ.dim),VQ,[Voxdir1out001,'IG_',LABTEMP,'_Pmask.nii']);
                    DynamicBC_write_NIFTI(reshape(THR4.*DQ,VQ.dim),VQ,[Voxdir1out001,'IG_',LABTEMP,'_Q_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR4).*DOR,VQ.dim),VQ,[Voxdir1out001,'IG_',LABTEMP,'_OR_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR4).*DORcihigh,VQ.dim),VQ,[Voxdir1out001,'IG_',LABTEMP,'_QRcihigh_masked.nii']);
                    DynamicBC_write_NIFTI(reshape(abs(THR4).*DORcilow,VQ.dim),VQ,[Voxdir1out001,'IG_',LABTEMP,'_ORcilow_masked.nii']);
                end
            end
        end
        
        if Para.Heat.StPermutation
            Voxdir2 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation',filesep];
            Voxdir2out05 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P05',filesep];
            mkdir(Voxdir2out05);
            Voxdir2out005 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P005',filesep];
            mkdir(Voxdir2out005);
            Voxdir2out01 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P01',filesep];
            mkdir(Voxdir2out01);
            Voxdir2out001 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P001',filesep];
            mkdir(Voxdir2out001);            
            files = dir([Voxdir2,'*.nii']);
            for ifiles = 1:length(files)
                [VP DP] = Dynamic_read_dir_NIFTI([Voxdir2,files(ifiles).name]);
                [pth nam ext] = fileparts([Voxdir2,files(ifiles).name]);
                
                THR1P = find(DP>0.95);
                THR1N = find(DP<0.05&DP>0);
                THR1 = zeros(size(DP));
                THR1(THR1P) = DP(THR1P);
                THR1(THR1N) = DP(THR1N)-1;
                DynamicBC_write_NIFTI(reshape(THR1,VQ.dim),VQ,[Voxdir2out05,nam,'marked_P.nii']);
                
                
                THR2P = find(DP>0.99);
                THR2N = find(DP<0.01&DP>0);
                THR2 = zeros(size(DP));
                THR2(THR2P) = DP(THR2P);
                THR2(THR2N) = DP(THR2N)-1;
                DynamicBC_write_NIFTI(reshape(THR2,VQ.dim),VQ,[Voxdir2out01,nam,'marked_P.nii']);
                
                
                THR3P = find(DP>0.995);
                THR3N = find(DP<0.005&DP>0);
                THR3 = zeros(size(DP));
                THR3(THR3P) = DP(THR3P);
                THR1(THR3N) = DP(THR3N)-1;
                DynamicBC_write_NIFTI(reshape(THR3,VQ.dim),VQ,[Voxdir2out005,nam,'marked_P.nii']);
                
                
                THR4P = find(DP>0.999);
                THR4N = find(DP<0.001&DP>0);
                THR4 = zeros(size(DP));
                THR4(THR4P) = DP(THR4P);
                THR4(THR4N) = DP(THR4N)-1;
                DynamicBC_write_NIFTI(reshape(THR4,VQ.dim),VQ,[Voxdir2out001,nam,'marked_P.nii']);
            end
            
            mats = dir([Voxdir2,'*.mat']);
            for imat = 1:length(mats)
                load([Voxdir2,mats(imat).name])
                % Raw
                PmapDir = [Voxdir2,files(imat).name];
                THRESHOLD = 0.05;
                C05 = ClusNum.ClusterNumB005;
                C05P = ClusNum.ClusterNumP005;
                C05N = ClusNum.ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,VQ,Voxdir2out05)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.ClusterNumB001;
                C01P = ClusNum.ClusterNumP001;
                C01N = ClusNum.ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,VQ,Voxdir2out01)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.ClusterNumB0005;
                C005P = ClusNum.ClusterNumP0005;
                C005N = ClusNum.ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,VQ,Voxdir2out005)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.ClusterNumB0001;
                C001P = ClusNum.ClusterNumP0001;
                C001N = ClusNum.ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,VQ,Voxdir2out001)
                
                % IG                
                PmapDir = [Voxdir2,files(imat+length(mats)).name];
                THRESHOLD = 0.05;
                C05 = ClusNum.IG_ClusterNumB005;
                C05P = ClusNum.IG_ClusterNumP005;
                C05N = ClusNum.IG_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,VQ,Voxdir2out05)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.IG_ClusterNumB001;
                C01P = ClusNum.IG_ClusterNumP001;
                C01N = ClusNum.IG_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,VQ,Voxdir2out01)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.IG_ClusterNumB0005;
                C005P = ClusNum.IG_ClusterNumP0005;
                C005N = ClusNum.IG_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,VQ,Voxdir2out005)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.IG_ClusterNumB0001;
                C001P = ClusNum.IG_ClusterNumP0001;
                C001N = ClusNum.IG_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,VQ,Voxdir2out001)
            end
        end
    end
    if Para.Heat.R
        if Para.Heat.StChisquare
            ROIdir = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise',filesep];
            ROIdirchi = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test',filesep];
            ROIdirchi05 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P05',filesep];
            ROIdirchi01 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P01',filesep];
            ROIdirchi005 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P005',filesep];
            ROIdirchi001 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P001',filesep];
            
            mkdir(ROIdirchi001);mkdir(ROIdirchi005);mkdir(ROIdirchi01);mkdir(ROIdirchi05);
            tempfold = dir(ROIdir);
            for itemp = 1:length(tempfold)-2       
                dirtemps = [ROIdir,tempfold(itemp+2).name,filesep];
                files1 = dir([dirtemps,'Group*markedAtlas.csv']);
                files2 = dir([dirtemps,'Group*markedBoth.csv']);
                files3 = dir([dirtemps,'Group*markedROI.csv']);
                if length(files1)>=2
                    C = nchoosek(1:length(files1),2);
                    dirtempschi = [ROIdirchi,tempfold(itemp+2).name,filesep];
                    dirtempschi05 = [ROIdirchi05,tempfold(itemp+2).name,filesep];
                    dirtempschi01 = [ROIdirchi01,tempfold(itemp+2).name,filesep];
                    dirtempschi005 = [ROIdirchi005,tempfold(itemp+2).name,filesep];
                    dirtempschi001 = [ROIdirchi001,tempfold(itemp+2).name,filesep];
                    mkdir(dirtempschi05);mkdir(dirtempschi01);mkdir(dirtempschi005);mkdir(dirtempschi001);
                    clear markatlas markboth markROI
                    for ifil = 1:length(files1)
                        markatlas(:,ifil) = csvread([dirtemps,files1(ifil).name]);
                        markboth(:,ifil) = csvread([dirtemps,files2(ifil).name]);
                        markROI(:,ifil) = csvread([dirtemps,files3(ifil).name]);
                    end
                    if length(files1)>2
                        P_anova = load([dirtempschi,'p_ANOVA_markedAtlas.csv']);
                        Q_anova = load([dirtempschi,'Q_ANOVA_markedAtlas.csv']);
                        P_anova(isnan(P_anova)) = 0;P_anova(isinf(P_anova)) = 0;
                        Q_anova(isnan(Q_anova)) = 0;Q_anova(isinf(Q_anova)) = 0;
                        P05 = P_anova<0.05&P_anova>0;
                        Q05 = Q_anova.*P05;
                        P01 = P_anova<0.01&P_anova>0;
                        Q01 = Q_anova.*P01;
                        P005 = P_anova<0.005&P_anova>0;
                        Q005 = Q_anova.*P005;
                        P001 = P_anova<0.001&P_anova>0;
                        Q001 = Q_anova.*P001;
                        csvwrite([dirtempschi05,'p_markedAtlas.csv'],P05)
                        csvwrite([dirtempschi05,'Q_markedAtlas.csv'],Q05)
                        csvwrite([dirtempschi01,'p_markedAtlas.csv'],P01)
                        csvwrite([dirtempschi01,'Q_markedAtlas.csv'],Q01)
                        csvwrite([dirtempschi005,'p_markedAtlas.csv'],P005)
                        csvwrite([dirtempschi005,'Q_markedAtlas.csv'],Q005)
                        csvwrite([dirtempschi001,'p_markedAtlas.csv'],P001)
                        csvwrite([dirtempschi001,'Q_markedAtlas.csv'],Q001)                        
                        
                        P_anova = load([dirtempschi,'p_ANOVA_markedBoth.csv']);
                        Q_anova = load([dirtempschi,'Q_ANOVA_markedBoth.csv']);
                        P_anova(isnan(P_anova)) = 0;P_anova(isinf(P_anova)) = 0;
                        Q_anova(isnan(Q_anova)) = 0;Q_anova(isinf(Q_anova)) = 0;
                        P05 = P_anova<0.05&P_anova>0;
                        Q05 = Q_anova.*P05;
                        P01 = P_anova<0.01&P_anova>0;
                        Q01 = Q_anova.*P01;
                        P005 = P_anova<0.005&P_anova>0;
                        Q005 = Q_anova.*P005;
                        P001 = P_anova<0.001&P_anova>0;
                        Q001 = Q_anova.*P001;
                        csvwrite([dirtempschi05,'p_markedBoth.csv'],P05)
                        csvwrite([dirtempschi05,'Q_markedBoth.csv'],Q05)
                        csvwrite([dirtempschi01,'p_markedBoth.csv'],P01)
                        csvwrite([dirtempschi01,'Q_markedBoth.csv'],Q01)
                        csvwrite([dirtempschi005,'p_markedBoth.csv'],P005)
                        csvwrite([dirtempschi005,'Q_markedBoth.csv'],Q005)
                        csvwrite([dirtempschi001,'p_markedBoth.csv'],P001)
                        csvwrite([dirtempschi001,'Q_markedBoth.csv'],Q001)
                                                
                        %%
                        P_anova = load([dirtempschi,'p_ANOVA_markedROI.csv']);
                        Q_anova = load([dirtempschi,'Q_ANOVA_markedROI.csv']);
                        P_anova(isnan(P_anova)) = 0;P_anova(isinf(P_anova)) = 0;
                        Q_anova(isnan(Q_anova)) = 0;Q_anova(isinf(Q_anova)) = 0;
                        P05 = P_anova<0.05&P_anova>0;
                        Q05 = Q_anova.*P05;
                        P01 = P_anova<0.01&P_anova>0;
                        Q01 = Q_anova.*P01;
                        P005 = P_anova<0.005&P_anova>0;
                        Q005 = Q_anova.*P005;
                        P001 = P_anova<0.001&P_anova>0;
                        Q001 = Q_anova.*P001;
                        csvwrite([dirtempschi05,'p_markedROI.csv'],P05)
                        csvwrite([dirtempschi05,'Q_markedROI.csv'],Q05)
                        csvwrite([dirtempschi01,'p_markedROI.csv'],P01)
                        csvwrite([dirtempschi01,'Q_markedROI.csv'],Q01)
                        csvwrite([dirtempschi005,'p_markedROI.csv'],P005)
                        csvwrite([dirtempschi005,'Q_markedROI.csv'],Q005)
                        csvwrite([dirtempschi001,'p_markedROI.csv'],P001)
                        csvwrite([dirtempschi001,'Q_markedROI.csv'],Q001)
                        
                    end
                    for j = 1:size(C,1)
                        LABname = ['Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2))];
                        
                        P = csvread([dirtempschi,'p_',LABname,'_markedAtlas.csv']);
                        Q = csvread([dirtempschi,'Q_',LABname,'_markedAtlas.csv']);
                        OR = csvread([dirtempschi,'OR_',LABname,'_markedAtlas.csv']);
                        ORCL = csvread([dirtempschi,'ORcilow_',LABname,'_markedAtlas.csv']);
                        ORCH = csvread([dirtempschi,'ORcihigh_',LABname,'_markedAtlas.csv']);
                        DIFdat = markatlas(:,C(j,1))-markatlas(:,C(j,2));
                        IND1P = P<0.05&P>0&DIFdat>0;
                        IND1N = P<0.05&P>0&DIFdat<0;
%                         IND1 = IND1P-IND1N;
                        IND1 = IND1P+IND1N;
                        csvwrite([dirtempschi05,'pmask_',LABname,'_markedAtlas.csv'],IND1.*P);   
                        csvwrite([dirtempschi05,'Q_marked_',LABname,'_markedAtlas.csv'],IND1.*Q); 
                        csvwrite([dirtempschi05,'OR_marked_',LABname,'_markedAtlas.csv'],abs(IND1).*OR);
                        csvwrite([dirtempschi05,'ORcilow_',LABname,'_markedAtlas.csv'],abs(IND1).*ORCL);
                        csvwrite([dirtempschi05,'ORcihigh_',LABname,'_markedAtlas.csv'],abs(IND1).*ORCH);
                        
                        IND2P = P<0.01&P>0&DIFdat>0;
                        IND2N = P<0.01&P>0&DIFdat<0;
%                         IND2 = IND2P-IND2N;
                        IND2 = IND2P+IND2N;
                        csvwrite([dirtempschi01,'pmask_',LABname,'_markedAtlas.csv'],IND2.*P);
                        csvwrite([dirtempschi01,'Q_marked_',LABname,'_markedAtlas.csv'],IND2.*Q);
                        csvwrite([dirtempschi01,'OR_marked_',LABname,'_markedAtlas.csv'],abs(IND2).*OR);
                        csvwrite([dirtempschi01,'ORcilow_',LABname,'_markedAtlas.csv'],abs(IND2).*ORCL);
                        csvwrite([dirtempschi01,'ORcihigh_',LABname,'_markedAtlas.csv'],abs(IND2).*ORCH);
                        IND3P = P<0.005&P>0&DIFdat>0;
                        IND3N = P<0.005&P>0&DIFdat<0;
%                         IND3 = IND3P-IND3N;
                        IND3 = IND3P+IND3N;
                        csvwrite([dirtempschi005,'pmask_',LABname,'_markedAtlas.csv'],IND3.*P);
                        csvwrite([dirtempschi005,'Q_marked_',LABname,'_markedAtlas.csv'],IND3.*Q);
                        csvwrite([dirtempschi005,'OR_marked_',LABname,'_markedAtlas.csv'],abs(IND3).*OR);
                        csvwrite([dirtempschi005,'ORcilow_',LABname,'_markedAtlas.csv'],abs(IND3).*ORCL);
                        csvwrite([dirtempschi005,'ORcihigh_',LABname,'_markedAtlas.csv'],abs(IND3).*ORCH);
                        IND4P = P<0.001&P>0&DIFdat>0;
                        IND4N = P<0.001&P>0&DIFdat<0;
%                         IND4 = IND4P-IND4N;
                        IND4 = IND4P+IND4N;
                        csvwrite([dirtempschi001,'pmask_',LABname,'_markedAtlas.csv'],IND4.*P);
                        csvwrite([dirtempschi001,'Q_marked_',LABname,'_markedAtlas.csv'],IND4.*Q);
                        csvwrite([dirtempschi001,'OR_marked_',LABname,'_markedAtlas.csv'],abs(IND4).*OR);
                        csvwrite([dirtempschi001,'ORcilow_',LABname,'_markedAtlas.csv'],abs(IND4).*ORCL);
                        csvwrite([dirtempschi001,'ORcihigh_',LABname,'_markedAtlas.csv'],abs(IND4).*ORCH);
                        %
                        
                        P = load([dirtempschi,'p_',LABname,'_markedBoth.csv']);
                        Q = load([dirtempschi,'Q_',LABname,'_markedBoth.csv']);
                        OR = load([dirtempschi,'OR_',LABname,'_markedBoth.csv']);
                        ORCL = load([dirtempschi,'ORcilow_',LABname,'_markedBoth.csv']);
                        ORCH = load([dirtempschi,'ORcihigh_',LABname,'_markedBoth.csv']);
                        DIFdat = markatlas(:,C(j,1))-markatlas(:,C(j,2));
                        IND1P = P<0.05&P>0&DIFdat>0;
                        IND1N = P<0.05&P>0&DIFdat<0;
%                         IND1 = IND1P-IND1N;
                        IND1 = IND1P+IND1N;
                        csvwrite([dirtempschi05,'pmask_',LABname,'_markedBoth.csv'],IND1.*P);
                        csvwrite([dirtempschi05,'Q_marked_',LABname,'_markedBoth.csv'],IND1.*Q);
                        csvwrite([dirtempschi05,'OR_marked_',LABname,'_markedBoth.csv'],abs(IND1).*OR);
                        csvwrite([dirtempschi05,'ORcilow_',LABname,'_markedBoth.csv'],abs(IND1).*ORCL);
                        csvwrite([dirtempschi05,'ORcihigh_',LABname,'_markedBoth.csv'],abs(IND1).*ORCH);
                        IND2P = P<0.01&P>0&DIFdat>0;
                        IND2N = P<0.01&P>0&DIFdat<0;
%                         IND2 = IND2P-IND2N;
                        IND2 = IND2P+IND2N;
                        csvwrite([dirtempschi01,'pmask_',LABname,'_markedBoth.csv'],IND2.*P);
                        csvwrite([dirtempschi01,'Q_marked_',LABname,'_markedBoth.csv'],IND2.*Q);
                        csvwrite([dirtempschi01,'OR_marked_',LABname,'_markedBoth.csv'],abs(IND2).*OR);
                        csvwrite([dirtempschi01,'ORcilow_',LABname,'_markedBoth.csv'],abs(IND2).*ORCL);
                        csvwrite([dirtempschi01,'ORcihigh_',LABname,'_markedBoth.csv'],abs(IND2).*ORCH);
                        IND3P = P<0.005&P>0&DIFdat>0;
                        IND3N = P<0.005&P>0&DIFdat<0;
%                         IND3 = IND3P-IND3N;
                        IND3 = IND3P+IND3N;
                        csvwrite([dirtempschi005,'pmask_',LABname,'_markedBoth.csv'],IND3.*P);
                        csvwrite([dirtempschi005,'Q_marked_',LABname,'_markedBoth.csv'],IND3.*Q);
                        csvwrite([dirtempschi005,'OR_marked_',LABname,'_markedBoth.csv'],abs(IND3).*OR);
                        csvwrite([dirtempschi005,'ORcilow_',LABname,'_markedBoth.csv'],abs(IND3).*ORCL);
                        csvwrite([dirtempschi005,'ORcihigh_',LABname,'_markedBoth.csv'],abs(IND3).*ORCH);
                        IND4P = P<0.001&P>0&DIFdat>0;
                        IND4N = P<0.001&P>0&DIFdat<0;
%                         IND4 = IND4P-IND4N;
                        IND4 = IND4P+IND4N;
                        csvwrite([dirtempschi001,'pmask_',LABname,'_markedBoth.csv'],IND4.*P);
                        csvwrite([dirtempschi001,'Q_marked_',LABname,'_markedBoth.csv'],IND4.*Q);
                        csvwrite([dirtempschi001,'OR_marked_',LABname,'_markedBoth.csv'],abs(IND4).*OR);
                        csvwrite([dirtempschi001,'ORcilow_',LABname,'_markedBoth.csv'],abs(IND4).*ORCL);
                        csvwrite([dirtempschi001,'ORcihigh_',LABname,'_markedBoth.csv'],abs(IND4).*ORCH);
                        %
                        
                        P = load([dirtempschi,'p_',LABname,'_markedROI.csv']);
                        Q = load([dirtempschi,'Q_',LABname,'_markedROI.csv']);
                        OR = load([dirtempschi,'OR_',LABname,'_markedROI.csv']);
                        ORCL = load([dirtempschi,'ORcilow_',LABname,'_markedROI.csv']);
                        ORCH = load([dirtempschi,'ORcihigh_',LABname,'_markedROI.csv']);
                        DIFdat = markatlas(:,C(j,1))-markatlas(:,C(j,2));
                        IND1P = P<0.05&P>0&DIFdat>0;
                        IND1N = P<0.05&P>0&DIFdat<0;
%                         IND1 = IND1P-IND1N;
                        IND1 = IND1P+IND1N;
                        csvwrite([dirtempschi05,'pmask_',LABname,'_markedROI.csv'],IND1.*P);
                        csvwrite([dirtempschi05,'Q_marked_',LABname,'_markedROI.csv'],IND1.*Q);
                        csvwrite([dirtempschi05,'OR_marked_',LABname,'_markedROI.csv'],abs(IND1).*OR);
                        csvwrite([dirtempschi05,'ORcilow_',LABname,'_markedROI.csv'],abs(IND1).*ORCL);
                        csvwrite([dirtempschi05,'ORcihigh_',LABname,'_markedROI.csv'],abs(IND1).*ORCH);
                        IND2P = P<0.01&P>0&DIFdat>0;
                        IND2N = P<0.01&P>0&DIFdat<0;
%                         IND2 = IND2P-IND2N;
                        IND2 = IND2P+IND2N;
                        csvwrite([dirtempschi01,'pmask_',LABname,'_markedROI.csv'],IND2.*P);
                        csvwrite([dirtempschi01,'Q_marked_',LABname,'_markedROI.csv'],IND2.*Q);
                        csvwrite([dirtempschi01,'OR_marked_',LABname,'_markedROI.csv'],abs(IND2).*OR);
                        csvwrite([dirtempschi01,'ORcilow_',LABname,'_markedROI.csv'],abs(IND2).*ORCL);
                        csvwrite([dirtempschi01,'ORcihigh_',LABname,'_markedROI.csv'],abs(IND2).*ORCH);
                        IND3P = P<0.005&P>0&DIFdat>0;
                        IND3N = P<0.005&P>0&DIFdat<0;
%                         IND3 = IND3P-IND3N;
                        IND3 = IND3P+IND3N;
                        csvwrite([dirtempschi005,'pmask_',LABname,'_markedROI.csv'],IND3.*P);
                        csvwrite([dirtempschi005,'Q_marked_',LABname,'_markedROI.csv'],IND3.*Q);
                        csvwrite([dirtempschi005,'OR_marked_',LABname,'_markedROI.csv'],abs(IND3).*OR);
                        csvwrite([dirtempschi005,'ORcilow_',LABname,'_markedROI.csv'],abs(IND3).*ORCL);
                        csvwrite([dirtempschi005,'ORcihigh_',LABname,'_markedROI.csv'],abs(IND3).*ORCH);
                        IND4P = P<0.001&P>0&DIFdat>0;
                        IND4N = P<0.001&P>0&DIFdat<0;
%                         IND4 = IND4P-IND4N;
                        IND4 = IND4P+IND4N;
                        csvwrite([dirtempschi001,'pmask_',LABname,'_markedROI.csv'],IND4.*P);
                        csvwrite([dirtempschi001,'Q_marked_',LABname,'_markedROI.csv'],IND4.*Q);
                        csvwrite([dirtempschi001,'OR_marked_',LABname,'_markedROI.csv'],abs(IND4).*OR);
                        csvwrite([dirtempschi001,'ORcilow_',LABname,'_markedROI.csv'],abs(IND4).*ORCL);
                        csvwrite([dirtempschi001,'ORcihigh_',LABname,'_markedROI.csv'],abs(IND4).*ORCH);
                    end
                end
                CSVfiletoNifti(dirtempschi001,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti(dirtempschi005,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti(dirtempschi01,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti(dirtempschi05,Vatl(1),Datl(:,itemp));
%                 dirtempschi001
%                 dirtempschi005
%                 dirtempschi01
%                 dirtempschi05
            end
        end
        if Para.Heat.StPermutation            
            ROIdirperm = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation',filesep];
            ROIdirperm05 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P05',filesep];
            ROIdirperm01 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P01',filesep];
            ROIdirperm005 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P005',filesep];
            ROIdirperm001 = [outdir,'HeatNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P001',filesep];
            mkdir(ROIdirperm001);mkdir(ROIdirperm005);mkdir(ROIdirperm01);mkdir(ROIdirperm05);
            tempfold = dir(ROIdirperm);
            for itemp = 1:length(tempfold)-2
                dirtemps = [ROIdirperm,tempfold(itemp+2).name,filesep];
                csvfile = dir([dirtemps,'*.csv']);
                mkdir([ROIdirperm05,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm01,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm005,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm001,tempfold(itemp+2).name,filesep]);
                for j = 1:length(csvfile)
                    PVAL = csvread([dirtemps,csvfile(j).name]);
                    Pindp = find(PVAL>0.95);
                    Pindn = find(PVAL<0.05&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm05,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.99);
                    Pindn = find(PVAL<0.01&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm01,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.995);
                    Pindn = find(PVAL<0.005&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm005,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.999);
                    Pindn = find(PVAL<0.001&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm001,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                end
            end
            
            CSVfiletoNifti([ROIdirperm001,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            CSVfiletoNifti([ROIdirperm005,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            CSVfiletoNifti([ROIdirperm01,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            CSVfiletoNifti([ROIdirperm05,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
        end
    end
    if Para.CentNum.R
        if Para.Heat.StChisquare
            ROIdir = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise',filesep];
            ROIdirchi = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test',filesep];
            ROIdirchi05 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P05',filesep];
            ROIdirchi01 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P01',filesep];
            ROIdirchi005 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P005',filesep];
            ROIdirchi001 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Chi2test_P001',filesep];
            mkdir(ROIdirchi001);mkdir(ROIdirchi005);mkdir(ROIdirchi01);mkdir(ROIdirchi05);
            
            tempfold = dir(ROIdir);
            for itemp = 1:length(tempfold)-2
                dirtemps = [ROIdir,tempfold(itemp+2).name,filesep];
                files1 = dir([dirtemps,'Group*.csv']);
                if length(files1)>=2
                    C = nchoosek(1:length(files1),2);
                    dirtempschi = [ROIdirchi,tempfold(itemp+2).name,filesep];
                    dirtempschi05 = [ROIdirchi05,tempfold(itemp+2).name,filesep];
                    dirtempschi01 = [ROIdirchi01,tempfold(itemp+2).name,filesep];
                    dirtempschi005 = [ROIdirchi005,tempfold(itemp+2).name,filesep];
                    dirtempschi001 = [ROIdirchi001,tempfold(itemp+2).name,filesep];
                    mkdir(dirtempschi05);mkdir(dirtempschi01);mkdir(dirtempschi005);mkdir(dirtempschi001);
                    clear markatlas markboth markROI
                    for ifil = 1:length(files1)
                        markatlas(:,ifil) = csvread([dirtemps,files1(ifil).name]);
                    end
                    if length(files1)>2
                        P_anova = csvread([dirtempschi,'p_ANOVA.csv']);
                        Q_anova = csvread([dirtempschi,'Q_ANOVA.csv']);
                        P_anova(isnan(P_anova)) = 0;P_anova(isinf(P_anova)) = 0;
                        Q_anova(isnan(Q_anova)) = 0;Q_anova(isinf(Q_anova)) = 0;
                        P05 = P_anova<0.05&P_anova>0;
                        Q05 = Q_anova.*P05;
                        P01 = P_anova<0.01&P_anova>0;
                        Q01 = Q_anova.*P01;
                        P005 = P_anova<0.005&P_anova>0;
                        Q005 = Q_anova.*P005;
                        P001 = P_anova<0.001&P_anova>0;
                        Q001 = Q_anova.*P001;
                        csvwrite([dirtempschi05,'p.csv'],P05)
                        csvwrite([dirtempschi05,'Q.csv'],Q05)
                        csvwrite([dirtempschi01,'p.csv'],P01)
                        csvwrite([dirtempschi01,'Q.csv'],Q01)
                        csvwrite([dirtempschi005,'p.csv'],P005)
                        csvwrite([dirtempschi005,'Q.csv'],Q005)
                        csvwrite([dirtempschi001,'p.csv'],P001)
                        csvwrite([dirtempschi001,'Q.csv'],Q001)
                    end
                    for j = 1:size(C,1)
                        LABname = ['Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2))];
                        
                        P = csvread([dirtempschi,'p_',LABname,'.csv']);
                        Q = csvread([dirtempschi,'Q_',LABname,'.csv']);
                        OR = csvread([dirtempschi,'OR_',LABname,'.csv']);
                        ORCL = csvread([dirtempschi,'ORcilow_',LABname,'.csv']);
                        ORCH = csvread([dirtempschi,'ORcihigh_',LABname,'.csv']);
                        DIFdat = markatlas(:,C(j,1))-markatlas(:,C(j,2));
                        IND1P = P<0.05&P>0&DIFdat>0;
                        IND1N = P<0.05&P>0&DIFdat<0;
%                         IND1 = IND1P-IND1N;
                        IND1 = IND1P+IND1N;
                        csvwrite([dirtempschi05,'pmask_',LABname,'.csv'],IND1.*P);
                        csvwrite([dirtempschi05,'Q_marked_',LABname,'.csv'],IND1.*Q);
                        csvwrite([dirtempschi05,'OR_marked_',LABname,'.csv'],abs(IND1).*OR);
                        csvwrite([dirtempschi05,'ORcilow_',LABname,'.csv'],abs(IND1).*ORCL);
                        csvwrite([dirtempschi05,'ORcihigh_',LABname,'.csv'],abs(IND1).*ORCH);
                        IND2P = P<0.01&P>0&DIFdat>0;
                        IND2N = P<0.01&P>0&DIFdat<0;
%                         IND2 = IND2P-IND2N;
                        IND2 = IND2P+IND2N;
                        csvwrite([dirtempschi01,'pmask_',LABname,'.csv'],IND2.*P);
                        csvwrite([dirtempschi01,'Q_marked_',LABname,'.csv'],IND2.*Q);
                        csvwrite([dirtempschi01,'OR_marked_',LABname,'.csv'],abs(IND2).*OR);
                        csvwrite([dirtempschi01,'ORcilow_',LABname,'.csv'],abs(IND2).*ORCL);
                        csvwrite([dirtempschi01,'ORcihigh_',LABname,'.csv'],abs(IND2).*ORCH);
                        IND3P = P<0.005&P>0&DIFdat>0;
                        IND3N = P<0.005&P>0&DIFdat<0;
%                         IND3 = IND3P-IND3N;
                        IND3 = IND3P+IND3N;
                        csvwrite([dirtempschi005,'pmask_',LABname,'.csv'],IND3.*P);
                        csvwrite([dirtempschi005,'Q_marked_',LABname,'.csv'],IND3.*Q);
                        csvwrite([dirtempschi005,'OR_marked_',LABname,'.csv'],abs(IND3).*OR);
                        csvwrite([dirtempschi005,'ORcilow_',LABname,'.csv'],abs(IND3).*ORCL);
                        csvwrite([dirtempschi005,'ORcihigh_',LABname,'.csv'],abs(IND3).*ORCH);
                        IND4P = P<0.001&P>0&DIFdat>0;
                        IND4N = P<0.001&P>0&DIFdat<0;
%                         IND4 = IND4P-IND4N;
                        IND4 = IND4P+IND4N;
                        csvwrite([dirtempschi001,'pmask_',LABname,'.csv'],IND4.*P);
                        csvwrite([dirtempschi001,'Q_marked_',LABname,'.csv'],IND4.*Q);
                        csvwrite([dirtempschi001,'OR_marked_',LABname,'.csv'],abs(IND4).*OR);
                        csvwrite([dirtempschi001,'ORcilow_',LABname,'.csv'],abs(IND4).*ORCL);
                        csvwrite([dirtempschi001,'ORcihigh_',LABname,'.csv'],abs(IND4).*ORCH);
                    end
                end
                
                CSVfiletoNifti(dirtempschi001,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti(dirtempschi005,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti(dirtempschi01,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti(dirtempschi05,Vatl(1),Datl(:,itemp));
            end
        end
        if Para.CentNum.StPermutation
            
            ROIdirperm = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation',filesep];
            ROIdirperm05 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P05',filesep];
            ROIdirperm01 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P01',filesep];
            ROIdirperm005 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P005',filesep];
            ROIdirperm001 = [outdir,'CenterNumber',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P001',filesep];
            mkdir(ROIdirperm001);mkdir(ROIdirperm005);mkdir(ROIdirperm01);mkdir(ROIdirperm05);
            tempfold = dir(ROIdirperm);
            for itemp = 1:length(tempfold)-2
                dirtemps = [ROIdirperm,tempfold(itemp+2).name,filesep];
                csvfile = dir([dirtemps,'*.csv']);
                mkdir([ROIdirperm05,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm01,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm005,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm001,tempfold(itemp+2).name,filesep]);
                for j = 1:length(csvfile)
                    PVAL = csvread([dirtemps,csvfile(j).name]);
                    Pindp = find(PVAL>0.95);
                    Pindn = find(PVAL<0.05&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm05,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.99);
                    Pindn = find(PVAL<0.01&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm01,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.995);
                    Pindn = find(PVAL<0.005&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm005,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.999);
                    Pindn = find(PVAL<0.001&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm001,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                end
            end
            
%             CSVfiletoNifti(ROIdirperm001,Vatl(1),Datl(:,itemp));
%             CSVfiletoNifti(ROIdirperm005,Vatl(1),Datl(:,itemp));
%             CSVfiletoNifti(ROIdirperm01,Vatl(1),Datl(:,itemp));
%             CSVfiletoNifti(ROIdirperm05,Vatl(1),Datl(:,itemp));
            
            CSVfiletoNifti([ROIdirperm001,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            CSVfiletoNifti([ROIdirperm005,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            CSVfiletoNifti([ROIdirperm01,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            CSVfiletoNifti([ROIdirperm05,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
        end
    end
    if Para.AffVolume.R
        if Para.AffVolume.StPermutation
            ROIdirperm = [outdir,'AffectedVolume',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation',filesep];
            ROIdirperm05 = [outdir,'AffectedVolume',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P05',filesep];
            ROIdirperm01 = [outdir,'AffectedVolume',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P01',filesep];
            ROIdirperm005 = [outdir,'AffectedVolume',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P005',filesep];
            ROIdirperm001 = [outdir,'AffectedVolume',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P001',filesep];
            mkdir(ROIdirperm001);mkdir(ROIdirperm005);mkdir(ROIdirperm01);mkdir(ROIdirperm05);
            tempfold = dir(ROIdirperm);
            for itemp = 1:length(tempfold)-2
                dirtemps = [ROIdirperm,tempfold(itemp+2).name,filesep];
                csvfile = dir([dirtemps,'*.csv']);
                mkdir([ROIdirperm05,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm01,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm005,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm001,tempfold(itemp+2).name,filesep]);
                for j = 1:length(csvfile)
                    PVAL = csvread([dirtemps,csvfile(j).name]);
                    Pindp = find(PVAL>0.95);
                    Pindn = find(PVAL<0.05&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm05,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.99);
                    Pindn = find(PVAL<0.01&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm01,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.995);
                    Pindn = find(PVAL<0.005&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm005,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.999);
                    Pindn = find(PVAL<0.001&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm001,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                end
                
%                 CSVfiletoNifti(ROIdirperm001,Vatl(1),Datl(:,itemp));
%                 CSVfiletoNifti(ROIdirperm005,Vatl(1),Datl(:,itemp));
%                 CSVfiletoNifti(ROIdirperm01,Vatl(1),Datl(:,itemp));
%                 CSVfiletoNifti(ROIdirperm05,Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm001,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm005,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm01,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm05,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            end
        end
    end
    if Para.WeiCent.V
        if Para.WeiCent.StPermutation
            Voxdir2 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation',filesep];
            Voxdir2out05 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P05',filesep];
            mkdir(Voxdir2out05);
            Voxdir2out005 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P005',filesep];
            mkdir(Voxdir2out005);
            Voxdir2out01 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P01',filesep];
            mkdir(Voxdir2out01);
            Voxdir2out001 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'VoxelWise_Permutation_P001',filesep];
            mkdir(Voxdir2out001);
            
            files = dir([Voxdir2,'*.nii']);
            for ifiles = 1:length(files)
                [VP DP] = Dynamic_read_dir_NIFTI([Voxdir2,files(ifiles).name]);
                [pth nam ext] = fileparts([Voxdir2,files(ifiles).name]);
                
                THR1P = find(DP>0.95);
                THR1N = find(DP<0.05&DP>0);
                THR1 = zeros(size(DP));
                THR1(THR1P) = DP(THR1P);
                THR1(THR1N) = DP(THR1N)-1;
                DynamicBC_write_NIFTI(reshape(THR1,VQ.dim),VQ,[Voxdir2out05,nam,'marked_P.nii']);
                
                
                THR2P = find(DP>0.99);
                THR2N = find(DP<0.01&DP>0);
                THR2 = zeros(size(DP));
                THR2(THR2P) = DP(THR2P);
                THR2(THR2N) = DP(THR2N)-1;
                DynamicBC_write_NIFTI(reshape(THR2,VQ.dim),VQ,[Voxdir2out01,nam,'marked_P.nii']);
                
                
                THR3P = find(DP>0.995);
                THR3N = find(DP<0.005&DP>0);
                THR3 = zeros(size(DP));
                THR3(THR3P) = DP(THR3P);
                THR1(THR3N) = DP(THR3N)-1;
                DynamicBC_write_NIFTI(reshape(THR3,VQ.dim),VQ,[Voxdir2out005,nam,'marked_P.nii']);
                
                
                THR4P = find(DP>0.999);
                THR4N = find(DP<0.001&DP>0);
                THR4 = zeros(size(DP));
                THR4(THR4P) = DP(THR4P);
                THR4(THR4N) = DP(THR4N)-1;
                DynamicBC_write_NIFTI(reshape(THR4,VQ.dim),VQ,[Voxdir2out001,nam,'marked_P.nii']);
            end
            
            mats = dir([Voxdir2,'*.mat']);
            
            for imat = 1:length(mats)                
                load([Voxdir2,mats(imat).name])
                % Raw
                PmapDir = [Voxdir2,files(imat).name];
                THRESHOLD = 0.05;
                C05 = ClusNum.ClusterNumB005;
                C05P = ClusNum.ClusterNumP005;
                C05N = ClusNum.ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,VQ,Voxdir2out05)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.ClusterNumB001;
                C01P = ClusNum.ClusterNumP001;
                C01N = ClusNum.ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,VQ,Voxdir2out01)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.ClusterNumB0005;
                C005P = ClusNum.ClusterNumP0005;
                C005N = ClusNum.ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,VQ,Voxdir2out005)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.ClusterNumB0001;
                C001P = ClusNum.ClusterNumP0001;
                C001N = ClusNum.ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,VQ,Voxdir2out001)
                %%
            end
            
        end
    end
    if Para.WeiCent.R
        if Para.WeiCent.StPermutation
            ROIdirperm = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation',filesep];
            ROIdirperm05 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P05',filesep];
            ROIdirperm01 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P01',filesep];
            ROIdirperm005 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P005',filesep];
            ROIdirperm001 = [outdir,'VolumeWeightedCenter',filesep,SetUpPara.ParaOut(i).LabName,filesep,'ROIWise_Permutation_P001',filesep];
            mkdir(ROIdirperm001);mkdir(ROIdirperm005);mkdir(ROIdirperm01);mkdir(ROIdirperm05);
            tempfold = dir(ROIdirperm);
            for itemp = 1:length(tempfold)-2
                dirtemps = [ROIdirperm,tempfold(itemp+2).name,filesep];
                csvfile = dir([dirtemps,'*.csv']);
                mkdir([ROIdirperm05,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm01,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm005,tempfold(itemp+2).name,filesep]);
                mkdir([ROIdirperm001,tempfold(itemp+2).name,filesep]);
                for j = 1:length(csvfile)
                    PVAL = csvread([dirtemps,csvfile(j).name]);
                    Pindp = find(PVAL>0.95);
                    Pindn = find(PVAL<0.05&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm05,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.99);
                    Pindn = find(PVAL<0.01&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm01,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.995);
                    Pindn = find(PVAL<0.005&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm005,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                    
                    Pindp = find(PVAL>0.999);
                    Pindn = find(PVAL<0.001&PVAL>0);
                    Pout = zeros(size(PVAL));
                    Pout(Pindp) = PVAL(Pindp);
                    Pout(Pindn) = PVAL(Pindn)-1;
                    csvwrite([ROIdirperm001,tempfold(itemp+2).name,filesep,csvfile(j).name],Pout)
                end
                
                CSVfiletoNifti([ROIdirperm001,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm005,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm01,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
                CSVfiletoNifti([ROIdirperm05,tempfold(itemp+2).name],Vatl(1),Datl(:,itemp));
            end
        end
    end    
    %%    
    if Para.WeiLoc.Vol
        if Para.WeiLoc.StPermutation
            outdirWeiVol = [outdir,'VolumeWeightedLocation',filesep,SetUpPara.ParaOut(i).LabName,filesep];
            RawDir = [outdirWeiVol,'VoxelWise_Permutation',filesep];
            files = dir([RawDir,'*.nii']);
            outdir1 = [outdirWeiVol,'VoxelWise_Permutation_P05',filesep];
            outdir2 = [outdirWeiVol,'VoxelWise_Permutation_P01',filesep];
            outdir3 = [outdirWeiVol,'VoxelWise_Permutation_P005',filesep];
            outdir4 = [outdirWeiVol,'VoxelWise_Permutation_P001',filesep];
            mkdir(outdir1);mkdir(outdir2);mkdir(outdir3);mkdir(outdir4);
            for j = 1:length(files)
                [v d] = Dynamic_read_dir_NIFTI([RawDir,files(j).name]);                
                [pth nam ext] = fileparts([RawDir,files(j).name]);
                
                THR1P = find(d>0.95);
                THR1N = find(d<0.05&d>0);
                THR1 = zeros(v.dim);
                THR1(THR1P) = d(THR1P);
                THR1(THR1N) = d(THR1N)-1;
                DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir1,nam,'marked_P.nii']);                
                THR1P = find(d>0.99);
                THR1N = find(d<0.01&d>0);
                THR1 = zeros(v.dim);
                THR1(THR1P) = d(THR1P);
                THR1(THR1N) = d(THR1N)-1;
                DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir2,nam,'marked_P.nii']);                
                THR1P = find(d>0.995);
                THR1N = find(d<0.005&d>0);
                THR1 = zeros(v.dim);
                THR1(THR1P) = d(THR1P);
                THR1(THR1N) = d(THR1N)-1;
                DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir3,nam,'marked_P.nii']);                
                THR1P = find(d>0.999);
                THR1N = find(d<0.001&d>0);
                THR1 = zeros(v.dim);
                THR1(THR1P) = d(THR1P);
                THR1(THR1N) = d(THR1N)-1;
                DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir4,nam,'marked_P.nii']);
            end
            
            mats = dir([RawDir,'*.mat']);
            
            for imat = 1:length(mats)                
                load([RawDir,mats(imat).name])
                extnam = mats(imat).name(1:end-11);
                % Raw
                PmapDir = [RawDir,extnam,'_mean_PermP.nii'];
                THRESHOLD = 0.05;
                C05 = ClusNum.mean_ClusterNumB005;
                C05P = ClusNum.mean_ClusterNumP005;
                C05N = ClusNum.mean_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.mean_ClusterNumB001;
                C01P = ClusNum.mean_ClusterNumP001;
                C01N = ClusNum.mean_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.mean_ClusterNumB0005;
                C005P = ClusNum.mean_ClusterNumP0005;
                C005N = ClusNum.mean_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.mean_ClusterNumB0001;
                C001P = ClusNum.mean_ClusterNumP0001;
                C001N = ClusNum.mean_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                
                %
                PmapDir = [RawDir,extnam,'_median_PermP.nii'];
                THRESHOLD = 0.05;
                C05 = ClusNum.media_ClusterNumB005;
                C05P = ClusNum.media_ClusterNumP005;
                C05N = ClusNum.media_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.media_ClusterNumB001;
                C01P = ClusNum.media_ClusterNumP001;
                C01N = ClusNum.media_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.media_ClusterNumB0005;
                C005P = ClusNum.media_ClusterNumP0005;
                C005N = ClusNum.media_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.media_ClusterNumB0001;
                C001P = ClusNum.media_ClusterNumP0001;
                C001N = ClusNum.media_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                
                %
                PmapDir = [RawDir,extnam,'_std_PermP.nii'];
                THRESHOLD = 0.05;
                C05 = ClusNum.std_ClusterNumB005;
                C05P = ClusNum.std_ClusterNumP005;
                C05N = ClusNum.std_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.std_ClusterNumB001;
                C01P = ClusNum.std_ClusterNumP001;
                C01N = ClusNum.std_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.std_ClusterNumB0005;
                C005P = ClusNum.std_ClusterNumP0005;
                C005N = ClusNum.std_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.std_ClusterNumB0001;
                C001P = ClusNum.std_ClusterNumP0001;
                C001N = ClusNum.std_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                %%               
                
                % IG
                PmapDir = [RawDir,'IG_',extnam,'_mean_PermP.nii'];
                THRESHOLD = 0.05;
                C05 = ClusNum.IG_mean_ClusterNumB005;
                C05P = ClusNum.IG_mean_ClusterNumP005;
                C05N = ClusNum.IG_mean_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.IG_mean_ClusterNumB001;
                C01P = ClusNum.IG_mean_ClusterNumP001;
                C01N = ClusNum.IG_mean_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.IG_mean_ClusterNumB0005;
                C005P = ClusNum.IG_mean_ClusterNumP0005;
                C005N = ClusNum.IG_mean_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.IG_mean_ClusterNumB0001;
                C001P = ClusNum.IG_mean_ClusterNumP0001;
                C001N = ClusNum.IG_mean_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                
                %
                PmapDir = [RawDir,'IG_',extnam,'_median_PermP.nii'];
                THRESHOLD = 0.05;
                C05 = ClusNum.IG_media_ClusterNumB005;
                C05P = ClusNum.IG_media_ClusterNumP005;
                C05N = ClusNum.IG_media_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.IG_media_ClusterNumB001;
                C01P = ClusNum.IG_media_ClusterNumP001;
                C01N = ClusNum.IG_media_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.IG_media_ClusterNumB0005;
                C005P = ClusNum.IG_media_ClusterNumP0005;
                C005N = ClusNum.IG_media_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.IG_media_ClusterNumB0001;
                C001P = ClusNum.IG_media_ClusterNumP0001;
                C001N = ClusNum.IG_media_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                
                %
                PmapDir = [RawDir,'IG_',extnam,'_std_PermP.nii'];
                THRESHOLD = 0.05;
                C05 = ClusNum.IG_std_ClusterNumB005;
                C05P = ClusNum.IG_std_ClusterNumP005;
                C05N = ClusNum.IG_std_ClusterNumN005;                
                CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                
                THRESHOLD = 0.01;
                C01 = ClusNum.IG_std_ClusterNumB001;
                C01P = ClusNum.IG_std_ClusterNumP001;
                C01N = ClusNum.IG_std_ClusterNumN001;                
                CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                
                THRESHOLD = 0.005;
                C005 = ClusNum.IG_std_ClusterNumB0005;
                C005P = ClusNum.IG_std_ClusterNumP0005;
                C005N = ClusNum.IG_std_ClusterNumN0005;                
                CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                
                THRESHOLD = 0.001;
                C001 = ClusNum.IG_std_ClusterNumB0001;
                C001P = ClusNum.IG_std_ClusterNumP0001;
                C001N = ClusNum.IG_std_ClusterNumN0001;                
                CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                %%
            end
        end
    end
    if Para.WeiLoc.Clin
        sfoldClin = dir([outdir,'ClincWeightedLocation*']);
        for ifold = 1:length(sfoldClin)
            if Para.WeiLoc.StPermutation
                outdirWeiVol = [outdir,sfoldClin(ifold).name,filesep,SetUpPara.ParaOut(i).LabName,filesep];
                RawDir = [outdirWeiVol,'VoxelWise_Permutation',filesep];
                files = dir([RawDir,'*.nii']);
                
                dir_PClus = [outdirWeiVol,'VoxelWise_PermutationCluster',filesep];
                PCluster_Group = dir([dir_PClus,'Group*']);
                permnums = load([dir_PClus,'permnums.mat']);
                
                outdir1 = [outdirWeiVol,'VoxelWise_Permutation_P05',filesep];
                outdir2 = [outdirWeiVol,'VoxelWise_Permutation_P01',filesep];
                outdir3 = [outdirWeiVol,'VoxelWise_Permutation_P005',filesep];
                outdir4 = [outdirWeiVol,'VoxelWise_Permutation_P001',filesep];
                mkdir(outdir1);mkdir(outdir2);mkdir(outdir3);mkdir(outdir4);
                for j = 1:length(files)
                    [v d] = Dynamic_read_dir_NIFTI([RawDir,files(j).name]);
                    [pth nam ext] = fileparts([RawDir,files(j).name]);
                    
                    THR1P = find(d>0.95);
                    THR1N = find(d<0.05&d>0);
                    THR1 = zeros(v.dim);
                    THR1(THR1P) = d(THR1P);
                    THR1(THR1N) = d(THR1N)-1;
                    DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir1,nam,'marked_P.nii']);
                    THR1P = find(d>0.99);
                    THR1N = find(d<0.01&d>0);
                    THR1 = zeros(v.dim);
                    THR1(THR1P) = d(THR1P);
                    THR1(THR1N) = d(THR1N)-1;
                    DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir2,nam,'marked_P.nii']);
                    THR1P = find(d>0.995);
                    THR1N = find(d<0.005&d>0);
                    THR1 = zeros(v.dim);
                    THR1(THR1P) = d(THR1P);
                    THR1(THR1N) = d(THR1N)-1;
                    DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir3,nam,'marked_P.nii']);
                    THR1P = find(d>0.999);
                    THR1N = find(d<0.001&d>0);
                    THR1 = zeros(v.dim);
                    THR1(THR1P) = d(THR1P);
                    THR1(THR1N) = d(THR1N)-1;
                    DynamicBC_write_NIFTI(reshape(THR1,v.dim),v,[outdir4,nam,'marked_P.nii']);
                end
                
                
                mats = dir([RawDir,'*.mat']);
                
                for imat = 1:length(mats)
                    load([RawDir,mats(imat).name])
                    extnam = mats(imat).name(1:end-11);
                    % Raw
                    PmapDir = [RawDir,extnam,'_mean_PermP.nii'];
                    THRESHOLD = 0.05;
                    C05 = ClusNum.mean_ClusterNumB005;
                    C05P = ClusNum.mean_ClusterNumP005;
                    C05N = ClusNum.mean_ClusterNumN005;
                    CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                    
                    THRESHOLD = 0.01;
                    C01 = ClusNum.mean_ClusterNumB001;
                    C01P = ClusNum.mean_ClusterNumP001;
                    C01N = ClusNum.mean_ClusterNumN001;
                    CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                    
                    THRESHOLD = 0.005;
                    C005 = ClusNum.mean_ClusterNumB0005;
                    C005P = ClusNum.mean_ClusterNumP0005;
                    C005N = ClusNum.mean_ClusterNumN0005;
                    CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                    
                    THRESHOLD = 0.001;
                    C001 = ClusNum.mean_ClusterNumB0001;
                    C001P = ClusNum.mean_ClusterNumP0001;
                    C001N = ClusNum.mean_ClusterNumN0001;
                    CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                    
                    %
                    PmapDir = [RawDir,extnam,'_median_PermP.nii'];
                    THRESHOLD = 0.05;
                    C05 = ClusNum.media_ClusterNumB005;
                    C05P = ClusNum.media_ClusterNumP005;
                    C05N = ClusNum.media_ClusterNumN005;
                    CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                    
                    THRESHOLD = 0.01;
                    C01 = ClusNum.media_ClusterNumB001;
                    C01P = ClusNum.media_ClusterNumP001;
                    C01N = ClusNum.media_ClusterNumN001;
                    CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                    
                    THRESHOLD = 0.005;
                    C005 = ClusNum.media_ClusterNumB0005;
                    C005P = ClusNum.media_ClusterNumP0005;
                    C005N = ClusNum.media_ClusterNumN0005;
                    CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                    
                    THRESHOLD = 0.001;
                    C001 = ClusNum.media_ClusterNumB0001;
                    C001P = ClusNum.media_ClusterNumP0001;
                    C001N = ClusNum.media_ClusterNumN0001;
                    CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                    
                    %
                    PmapDir = [RawDir,extnam,'_std_PermP.nii'];
                    THRESHOLD = 0.05;
                    C05 = ClusNum.std_ClusterNumB005;
                    C05P = ClusNum.std_ClusterNumP005;
                    C05N = ClusNum.std_ClusterNumN005;
                    CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                    
                    THRESHOLD = 0.01;
                    C01 = ClusNum.std_ClusterNumB001;
                    C01P = ClusNum.std_ClusterNumP001;
                    C01N = ClusNum.std_ClusterNumN001;
                    CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                    
                    THRESHOLD = 0.005;
                    C005 = ClusNum.std_ClusterNumB0005;
                    C005P = ClusNum.std_ClusterNumP0005;
                    C005N = ClusNum.std_ClusterNumN0005;
                    CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                    
                    THRESHOLD = 0.001;
                    C001 = ClusNum.std_ClusterNumB0001;
                    C001P = ClusNum.std_ClusterNumP0001;
                    C001N = ClusNum.std_ClusterNumN0001;
                    CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                    %%
                    
                    
                    % IG
                    PmapDir = [RawDir,'IG_',extnam,'_mean_PermP.nii'];
                    THRESHOLD = 0.05;
                    C05 = ClusNum.IG_mean_ClusterNumB005;
                    C05P = ClusNum.IG_mean_ClusterNumP005;
                    C05N = ClusNum.IG_mean_ClusterNumN005;
                    CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                    
                    THRESHOLD = 0.01;
                    C01 = ClusNum.IG_mean_ClusterNumB001;
                    C01P = ClusNum.IG_mean_ClusterNumP001;
                    C01N = ClusNum.IG_mean_ClusterNumN001;
                    CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                    
                    THRESHOLD = 0.005;
                    C005 = ClusNum.IG_mean_ClusterNumB0005;
                    C005P = ClusNum.IG_mean_ClusterNumP0005;
                    C005N = ClusNum.IG_mean_ClusterNumN0005;
                    CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                    
                    THRESHOLD = 0.001;
                    C001 = ClusNum.IG_mean_ClusterNumB0001;
                    C001P = ClusNum.IG_mean_ClusterNumP0001;
                    C001N = ClusNum.IG_mean_ClusterNumN0001;
                    CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                    
                    %
                    PmapDir = [RawDir,'IG_',extnam,'_median_PermP.nii'];
                    THRESHOLD = 0.05;
                    C05 = ClusNum.IG_media_ClusterNumB005;
                    C05P = ClusNum.IG_media_ClusterNumP005;
                    C05N = ClusNum.IG_media_ClusterNumN005;
                    CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                    
                    THRESHOLD = 0.01;
                    C01 = ClusNum.IG_media_ClusterNumB001;
                    C01P = ClusNum.IG_media_ClusterNumP001;
                    C01N = ClusNum.IG_media_ClusterNumN001;
                    CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                    
                    THRESHOLD = 0.005;
                    C005 = ClusNum.IG_media_ClusterNumB0005;
                    C005P = ClusNum.IG_media_ClusterNumP0005;
                    C005N = ClusNum.IG_media_ClusterNumN0005;
                    CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                    
                    THRESHOLD = 0.001;
                    C001 = ClusNum.IG_media_ClusterNumB0001;
                    C001P = ClusNum.IG_media_ClusterNumP0001;
                    C001N = ClusNum.IG_media_ClusterNumN0001;
                    CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                    
                    %
                    PmapDir = [RawDir,'IG_',extnam,'_std_PermP.nii'];
                    THRESHOLD = 0.05;
                    C05 = ClusNum.IG_std_ClusterNumB005;
                    C05P = ClusNum.IG_std_ClusterNumP005;
                    C05N = ClusNum.IG_std_ClusterNumN005;
                    CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,v,outdir1)
                    
                    THRESHOLD = 0.01;
                    C01 = ClusNum.IG_std_ClusterNumB001;
                    C01P = ClusNum.IG_std_ClusterNumP001;
                    C01N = ClusNum.IG_std_ClusterNumN001;
                    CFA_MakeRes_Clustercorrected(C01,C01P,C01N,PmapDir,THRESHOLD,v,outdir2)
                    
                    THRESHOLD = 0.005;
                    C005 = ClusNum.IG_std_ClusterNumB0005;
                    C005P = ClusNum.IG_std_ClusterNumP0005;
                    C005N = ClusNum.IG_std_ClusterNumN0005;
                    CFA_MakeRes_Clustercorrected(C005,C005P,C005N,PmapDir,THRESHOLD,v,outdir3)
                    
                    THRESHOLD = 0.001;
                    C001 = ClusNum.IG_std_ClusterNumB0001;
                    C001P = ClusNum.IG_std_ClusterNumP0001;
                    C001N = ClusNum.IG_std_ClusterNumN0001;
                    CFA_MakeRes_Clustercorrected(C001,C001P,C001N,PmapDir,THRESHOLD,v,outdir4)
                    %%
                end
            end
        end
    end
    %%
end
if Para.VLSM.ana
    outdirVLSM = [outdir,'VLSM',filesep];
    sfoldvlsm = dir(outdirVLSM);
    for j = 1:length(sfoldvlsm)-2
        dirtemp = [outdirVLSM,sfoldvlsm(j+2).name,filesep,'VoxelWise',filesep];
        
        outdirtemp1 = [outdirVLSM,sfoldvlsm(j+2).name,filesep,'VoxelWise_p05',filesep];
        outdirtemp2 = [outdirVLSM,sfoldvlsm(j+2).name,filesep,'VoxelWise_p01',filesep];
        outdirtemp3 = [outdirVLSM,sfoldvlsm(j+2).name,filesep,'VoxelWise_p005',filesep];
        outdirtemp4 = [outdirVLSM,sfoldvlsm(j+2).name,filesep,'VoxelWise_p001',filesep];
        mkdir(outdirtemp1);mkdir(outdirtemp2);mkdir(outdirtemp3);mkdir(outdirtemp4);
        [vt,dt] = Dynamic_read_dir_NIFTI([dirtemp,'Tmap.nii']);
        [vp,dp] = Dynamic_read_dir_NIFTI([dirtemp,'Pmap.nii']);
        [vz,dz] = Dynamic_read_dir_NIFTI([dirtemp,'Zmap.nii']);
        dt1 = dt.*(dp<0.05);dt1(isnan(dt1)) = 0;dt1(isinf(dt1)) = 0;
        dt2 = dt.*(dp<0.01);dt2(isnan(dt2)) = 0;dt2(isinf(dt2)) = 0;
        dt3 = dt.*(dp<0.005);dt3(isnan(dt3)) = 0;dt3(isinf(dt3)) = 0;
        dt4 = dt.*(dp<0.001);dt4(isnan(dt4)) = 0;dt4(isinf(dt4)) = 0;
        DynamicBC_write_NIFTI(reshape(dt1,vt.dim),vt,[outdirtemp1,'Tmap_mask.nii']);
        DynamicBC_write_NIFTI(reshape(dt2,vt.dim),vt,[outdirtemp2,'Tmap_mask.nii']);
        DynamicBC_write_NIFTI(reshape(dt3,vt.dim),vt,[outdirtemp3,'Tmap_mask.nii']);
        DynamicBC_write_NIFTI(reshape(dt4,vt.dim),vt,[outdirtemp4,'Tmap_mask.nii']);
        
        dz1 = dz.*(dp<0.05);dz1(isnan(dz1)) = 0;dz1(isinf(dz1)) = 0;
        dz2 = dz.*(dp<0.01);dz2(isnan(dz2)) = 0;dz2(isinf(dz2)) = 0;
        dz3 = dz.*(dp<0.005);dz3(isnan(dz3)) = 0;dz3(isinf(dz3)) = 0;
        dz4 = dz.*(dp<0.001);dz4(isnan(dz4)) = 0;dz4(isinf(dz4)) = 0;
        DynamicBC_write_NIFTI(reshape(dz1,vz.dim),vt,[outdirtemp1,'Zmap_mask.nii']);
        DynamicBC_write_NIFTI(reshape(dz2,vz.dim),vt,[outdirtemp2,'Zmap_mask.nii']);
        DynamicBC_write_NIFTI(reshape(dz3,vz.dim),vt,[outdirtemp3,'Zmap_mask.nii']);
        DynamicBC_write_NIFTI(reshape(dz4,vz.dim),vt,[outdirtemp4,'Zmap_mask.nii']);
        
        ClusNum = VLSM_clusterPermutation(dirtemp);
        save([outdirVLSM,sfoldvlsm(j+2).name,filesep,'ClusNum.mat'],'ClusNum')
%         load([outdirVLSM,sfoldvlsm(j+2).name,filesep,'ClusNum.mat'],'ClusNum')
        PmapDir = [dirtemp,'Pmap.nii'];
        TmapDir = [dirtemp,'Tmap.nii'];
        ZmapDir = [dirtemp,'Zmap.nii'];
        
        THRESHOLD = 0.05;
        C05 = ClusNum.ClusterNumB005;
        C05P = ClusNum.ClusterNumP005;
        C05N = ClusNum.ClusterNumN005;
        CFA_MakeRes_Clustercorrected_vlsm(C05,C05P,C05N,PmapDir,TmapDir,ZmapDir,THRESHOLD,vt,outdirtemp1)
        
        
        THRESHOLD = 0.01;
        C01 = ClusNum.ClusterNumB001;
        C01P = ClusNum.ClusterNumP001;
        C01N = ClusNum.ClusterNumN001;
        CFA_MakeRes_Clustercorrected_vlsm(C01,C01P,C01N,PmapDir,TmapDir,ZmapDir,THRESHOLD,vt,outdirtemp2)
        
        
        THRESHOLD = 0.005;
        C005 = ClusNum.ClusterNumB0005;
        C005P = ClusNum.ClusterNumP0005;
        C005N = ClusNum.ClusterNumN0005;
        CFA_MakeRes_Clustercorrected_vlsm(C005,C005P,C005N,PmapDir,TmapDir,ZmapDir,THRESHOLD,vt,outdirtemp3)
        
        
        THRESHOLD = 0.001;
        C001 = ClusNum.ClusterNumB0001;
        C001P = ClusNum.ClusterNumP0001;
        C001N = ClusNum.ClusterNumN0001;
        CFA_MakeRes_Clustercorrected_vlsm(C001,C001P,C001N,PmapDir,TmapDir,ZmapDir,THRESHOLD,vt,outdirtemp4)
        
    end
end
end