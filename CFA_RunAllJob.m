function CFA_RunAllJob(SetUpPara,Template,Para)
Indir = SetUpPara.Indir;
Outdir = SetUpPara.Outdir;
filp = which('CFA_RunAllJob.m');
[pth,nam,ext] = fileparts(filp);
spmjobdir = [pth,filesep,'NormJobmat',filesep];
Atlasdir = [pth,filesep,'TemplateUsed',filesep];
voxelsize = SetUpPara.OutVsize;
BrainMaskdir = [pth,filesep,'brainmask.nii'];

IGper = Para.IgP;

%% Normalize
if SetUpPara.NormOpt
    outdir = [Outdir,filesep,'Normalized',filesep];
    mkdir(outdir);
    sfold = dir(Indir);
    if SetUpPara.Normtype1 % epimode
        for i = 1:length(sfold)-2
            files = dir([Indir,filesep,sfold(i+2).name,filesep,'*.nii']);
%             inum = 1;
            if length(files)>2
                for j = 1:length(files)
                    if ~strcmp(files(j).name,'image.nii')&&~strcmp(files(j).name,'ROI.nii')
                        delete([Indir,filesep,sfold(i+2).name,filesep,files(j).name]);
                    end
                end
            end
        end
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            ClincFreq_Norm2MNI_EPI(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
%         CFA_CheckNormROI_batch(outdir,0);
    elseif SetUpPara.Normtype2 % T1mode
        
        for i = 1:length(sfold)-2
            files = dir([Indir,filesep,sfold(i+2).name,filesep,'*.nii']);
%             inum = 1;
            if length(files)>2
                for j = 1:length(files)
                    if ~strcmp(files(j).name,'image.nii')&&~strcmp(files(j).name,'ROI.nii')
                        delete([Indir,filesep,sfold(i+2).name,filesep,files(j).name]);
                    end
                end
            end
        end
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            %             ClincFreq_Norm2MNI_T1(indat,inroi,outdir,outname,voxelsize,spmjobdir)
            ClincFreq_Norm2MNI_T1_clinc(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
%         CFA_CheckNormROI_batch(outdir,0);
    elseif SetUpPara.Normtype3 % with T1mode
        
        for i = 1:length(sfold)-2
            files = dir([Indir,filesep,sfold(i+2).name,filesep,'*.nii']);
%             inum = 1;
            if length(files)>2
                for j = 1:length(files)
                    if ~strcmp(files(j).name,'T1.nii')&&~strcmp(files(j).name,'ROI.nii')&&~strcmp(files(j).name,'OtherMode.nii')
                        delete([Indir,filesep,sfold(i+2).name,filesep,files(j).name]);
                    end
                end
            end
        end
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'OtherMode.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            indatT1 = [Indir,filesep,sfold(i+2).name,filesep,'T1.nii'];
            outname = sfold(i+2).name;
            %             ClincFreq_Norm2MNI_withT1(indat,inroi,indatT1,outdir,outname,voxelsize,spmjobdir)
            ClincFreq_Norm2MNI_withT1_clinc(indat,inroi,indatT1,outdir,outname,voxelsize,spmjobdir)
        end        
%         CFA_CheckNormROI_batch(outdir,0);
    elseif SetUpPara.Normtype4 % T2mode
        
        for i = 1:length(sfold)-2
            files = dir([Indir,filesep,sfold(i+2).name,filesep,'*.nii']);
%             inum = 1;
            if length(files)>2
                for j = 1:length(files)
                    if ~strcmp(files(j).name,'image.nii')&&~strcmp(files(j).name,'ROI.nii')
                        delete([Indir,filesep,sfold(i+2).name,filesep,files(j).name]);
                    end
                end
            end
        end
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            ClincFreq_Norm2MNI_T2(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
%         CFA_CheckNormROI_batch(outdir,0);
    elseif SetUpPara.Normtype5 % CTmode
        
        for i = 1:length(sfold)-2
            files = dir([Indir,filesep,sfold(i+2).name,filesep,'*.nii']);
%             inum = 1;
            if length(files)>2
                for j = 1:length(files)
                    if ~strcmp(files(j).name,'image.nii')&&~strcmp(files(j).name,'ROI.nii')
                        delete([Indir,filesep,sfold(i+2).name,filesep,files(j).name]);
                    end
                end
            end
        end
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            %             ClincFreq_Norm2MNI_CT(indat,inroi,outdir,outname,voxelsize,spmjobdir)
            ClincFreq_Norm2MNI_CT_clinc(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
%         CFA_CheckNormROI_batch(outdir,1);
    end
    
    INDIR = [Outdir,filesep,'NormROI',filesep];
    mkdir(INDIR);
    copyfile([outdir,'ROI*.nii'],INDIR);
    INDIRBG = [Outdir,filesep,'NormBG',filesep];
    mkdir(INDIRBG);
    copyfile([outdir,'image*.nii'],INDIRBG);
    if ~isempty(dir([outdir,'T1*.nii']))
        INDIRT1bg = [Outdir,filesep,'NormBGT1',filesep];
        mkdir(INDIRT1bg);
        copyfile([outdir,'T1*.nii'],INDIRT1bg);
    end
    
%     Indir
%     Outdir
    filesROI = dir([outdir,'ROI*.nii']);
    if length(filesROI)~=(length(sfold)-2)
        for i = 1:length(sfold)-2
            sfoldlist{i,1} = i;
            sfoldlist{i,2} = sfold(i+2).name;
        end
%         LABexistind = 1;
        for i = 1:length(filesROI)
            nametemp = filesROI(i).name;
            inds = find(nametemp=='_');
            indslen = length(inds);
            NameTemp = nametemp(inds(1)+1:inds(indslen-1)-1);
%             elab = 0;
            for j = 1:length(sfold)-2
                if strcmp(NameTemp,sfoldlist{j,2})
                    LABexistind(i,1) = sfoldlist{j,1};
                    break
                end
            end            
        end
        SetUpParaOld = SetUpPara;
        for i = 2:length(SetUpPara.ParaOut)
            DATALAB = SetUpPara.ParaOut(i).datalab;
            for j = 1:length(DATALAB)
                datalabind = DATALAB{j};
                datasep = intersect(LABexistind,datalabind);
                DATALABnew{j} = datasep;
            end
            SetUpPara.ParaOut(i).datalab = DATALABnew;
        end
        save([SetUpPara.Outdir,filesep,'SetUpParaOld.mat'],'SetUpParaOld');
        save([SetUpPara.Outdir,filesep,'SetUpParaNew.mat'],'SetUpPara');
        save([SetUpPara.Outdir,filesep,'ExistInd.mat'],'LABexistind');
        
    end
    if ~SetUpPara.Normtype5
        CFA_CheckNormROI_batch(Outdir,0)
    else SetUpPara.Normtype5
        CFA_CheckNormROI_batch(Outdir,1) 
    end
else
    INDIR = [Indir,filesep];
    if Para.Other.Hist
        Indirtemp = uigetdir(pwd,'Histgram from data');
        INDIRBG = [Indirtemp,filesep];
    end
end

% dynamicBC_Reslice(AtlasTemp{i,1},AtlasTempU{i,1},voxelsize,0,[INDIR,filesep,ROIfil(1).name]);
Rfile = dir([INDIR,filesep,'*.nii']);
dynamicBC_Reslice(BrainMaskdir,[Outdir,filesep,'Bmask.nii'],voxelsize*[1 1 1],0,[INDIR,filesep,Rfile(1).name]);
[vbmask dbmask0] = Dynamic_read_dir_NIFTI([Outdir,filesep,'Bmask.nii']);
dbmask = reshape(dbmask0,vbmask.dim);

%
SetUpPara.ParaOut(1).datalab{1} = 1:length(Rfile);
%%
atlasnum = 0;
if Template.Fivetissue
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_FiveTissue.nii'];
    load([Atlasdir,'AtlasU_FiveTissue.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'FiveTissue';
end
if Template.FivetissueLR
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_FiveTissueLR.nii'];
    load([Atlasdir,'AtlasU_FiveTissueLR.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'FiveTissueLR';
end
if Template.JHUtotal
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_JHUtotal.nii'];
    load([Atlasdir,'AtlasU_JHUtotal.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'JHUtotal';
end
if Template.LPBA40
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_LPBA40.nii'];
    load([Atlasdir,'AtlasU_LPBA40.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'LPBA40';
end
if Template.AAL
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_AAL.nii'];
    load([Atlasdir,'AtlasU_AAL.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'AAL';
end
if Template.HOAcortex
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_HOAcortex.nii'];
    load([Atlasdir,'AtlasU_HOAcortex.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'HOAcortex';
end
if Template.HOAsubcortex
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_HOAsubcortex.nii'];
    load([Atlasdir,'AtlasU_HOAsubcortex.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'HOAsubcortex';
end
if Template.MesulamHOA
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_MesulamHOA.nii'];
    load([Atlasdir,'AtlasU_MesulamHOA.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'MesulamHOA';
end
if Template.JHUwhite
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_JHUwhite.nii'];
    load([Atlasdir,'AtlasU_JHUwhite.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'JHUwhite';
end
if Template.ICBMwhite
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_ICBMwhite.nii'];
    load([Atlasdir,'AtlasU_ICBMwhite.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'ICBMwhite';
end
if Template.Ventricle
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_Ventricle.nii'];
    load([Atlasdir,'AtlasU_Ventricle.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'Ventricle';
end
if Template.Arterial
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_Arterial.nii'];
    load([Atlasdir,'AtlasU_Arterial.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'Arteiral';
end
if Template.Other
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = Template.OtherROI;
    load(Template.OtherROIname);
    AtlasTemp{atlasnum,2} = ROItab;
    [pths,nams,exts] = fileparts(Template.OtherROI);
    AtlasTemp{atlasnum,3} = nams;
end

%%
if ~isempty(whos('AtlasTemp'))
    OutdirAtlas = [Outdir,filesep,'ReslicedAtlas',filesep];
    if ~isfolder(OutdirAtlas)
        mkdir(OutdirAtlas)
    end
    ROIfil = dir([INDIR,filesep,'*.nii']);
    for i = 1:size(AtlasTemp,1)
        tempdir = AtlasTemp{i,1};
        [pt,na,ex] = fileparts(tempdir);
        if Template.Other&&i==size(AtlasTemp,1)
            AtlasTempU{i,1} = [OutdirAtlas,filesep,'AtlasU_',na,ex];
        else
            AtlasTempU{i,1} = [OutdirAtlas,filesep,na,ex];
        end
        AtlasTempU{i,2} = AtlasTemp{i,2};
        AtlasTempU{i,3} = AtlasTemp{i,3};
        dynamicBC_Reslice(AtlasTemp{i,1},AtlasTempU{i,1},voxelsize,0,[INDIR,filesep,ROIfil(1).name]);
    end
    
    for i = 1:size(AtlasTempU,1)
        [VT DT] = Dynamic_read_dir_NIFTI(AtlasTempU{i,1});
        DT = round(DT);
        ATMark = AtlasTempU{i,2};
        for j = 1:size(ATMark,1)
            labval(j,1) = ATMark{j,1};
            labname{j,1} = ATMark{j,2};
        end
        indTab = find(labval>0);
        AtlasTempU{i,4} = labval(indTab);
        AtlasTempU{i,5} = labname(indTab);
        for j = 1:length(indTab)
            AtlasTempInd{j,1} = find(DT==labval(indTab(j)));
        end
        AtlasTempU{i,6} = AtlasTempInd;
        clear labval labname AtlasTempInd
    end
    if atlasnum==0
        AtlasTempU = [];
    end
end


if ~Template.Ventricle&&Para.Other.D2V
    OutdirAtlas = [Outdir,filesep,'ReslicedAtlas',filesep];
    if ~isfolder(OutdirAtlas)
        mkdir(OutdirAtlas)
    end
    ROIfil = dir([INDIR,filesep,'*.nii']);
    dynamicBC_Reslice([Atlasdir,'AtlasU_Ventricle.nii'],...
        [Outdir,filesep,'ReslicedAtlas',filesep,'AtlasU_Ventricle.nii'],voxelsize,...
        0,[INDIR,filesep,ROIfil(1).name]);
    copyfile([Atlasdir,'AtlasU_Ventricle.mat'],[Outdir,filesep,'ReslicedAtlas',filesep,'AtlasU_Ventricle.mat']);
elseif Para.Other.D2V
    copyfile([Atlasdir,'AtlasU_Ventricle.mat'],[Outdir,filesep,'ReslicedAtlas',filesep,'AtlasU_Ventricle.mat']);    
end

%% Load informations
if Para.WeiLoc.Clin 
    answ = inputdlg('Clinc Variable Number','WeiLocClinc',1,{'1'});
    ClincNum = str2num(answ{1});
    clear answ
    for icn = 1:ClincNum
        answ = inputdlg(['Var',num2str(icn),' Name:'],['Var',num2str(icn),' Name:'],1);
        ClincWei(icn).Name = answ{1};
        clear answ;
        [pg,pth] = uigetfile('*.txt',['Please Select Var ',ClincWei(icn).Name]);
        vartemp = load(fullfile(pth,pg));
        if exist('LABexistind')
            ClincWei(icn).Var = vartemp(LABexistind);
        else
            ClincWei(icn).Var = vartemp;
        end
    end
end
if Para.VLSM.ana
    ParaVLSM = {};
    ParaNum = 1;
    if Para.VLSM.VarCon
        for i = 2:length(SetUpPara.ParaOut)
            if SetUpPara.ParaOut(i).dattype~=1
                ParaVLSM{ParaNum,1} = SetUpPara.ParaOut(i).LabName;
                pvlsm = load(SetUpPara.ParaOut(i).Labtxt);
                if exist('LABexistind')
                    ParaVLSM{ParaNum,2} = pvlsm(LABexistind);
                else
                    ParaVLSM{ParaNum,2} = pvlsm;
                end
                ParaNum = ParaNum+1;
            end
        end
    end
    if Para.VLSM.VarNew
        valnum00 = inputdlg('tablelistnum','VLSMvar',1,{'1'});
        valnum = str2num(valnum00{1});
        for i = 1:valnum
            varnam00 = inputdlg(['Var ',num2str(i),' Name'],['Var ',num2str(i),' Name'],1);
            ParaVLSM{ParaNum,1} = varnam00{1};
            [fil,pth] = uigetfile('*.txt',[varnam00{1},' selection']);
            pvlsm = load(fullfile(pth,fil));
            if exist('LABexistind')
                ParaVLSM{ParaNum,2} = pvlsm(LABexistind);
            else
                ParaVLSM{ParaNum,2} = pvlsm;
            end
%             ParaVLSM{ParaNum,2} = load(fullfile(pth,fil));
            ParaNum = ParaNum+1;
        end
    end
    if Para.VLSM.COV     
        COVNum = 1;
        covnum00 = inputdlg('VLSM：number of var of no interest','No.ofVarnoInterest',1,{'1'});
        covnum = str2num(covnum00{1});
        for i = 1:covnum
            [fil,pth] = uigetfile('*.txt',['VLSM： cov ',num2str(i),' selection']);
            COVVLSM(:,COVNum) = load(fullfile(pth,fil));
            COVNum = COVNum+1;
        end
        if exist('LABexistind')
            COVVLSMO = COVVLSM;
            clear COVVLSM;
            COVVLSM = COVVLSMO(LABexistind,:);
        end
    else
        COVVLSM = [];
    end
end

%%
if Para.VLSM.ana
    
    OutVLSM = [Outdir,filesep,'VLSM',filesep];
    mkdir(OutVLSM);
    disp('Now Calculate VLSM!');
    ClincFreq_VLSM(INDIR,OutVLSM,SetUpPara,Para,dbmask,ParaVLSM,COVVLSM);
    disp('VLSM Calculatation finished!')
end

%%
if Para.Heat.V||Para.Heat.R
    OutHeatdir = [Outdir,filesep,'HeatNumber',filesep];
    mkdir(OutHeatdir);
    disp('Now Calculate HeatNumber!')
    ClincFreq_CalHeatNum(INDIR,OutHeatdir,SetUpPara,Para,AtlasTempU,dbmask)
    % stat needed
    ClincFreq_HeatNumStat(INDIR,OutHeatdir,SetUpPara,Para,AtlasTempU,dbmask);
    disp('HeatNumber Calculation finished')
end
%%
%%
if Para.CentNum.R||Para.WeiCent.V||Para.WeiCent.R
    outdirCent = [Outdir,filesep,'ForCentRelated',filesep];
    if isempty(dir(outdirCent))
        mkdir(outdirCent)
        ClincFreq_CentPointForVolume_pre(INDIR,outdirCent,Para,spmjobdir);
    else
        if Para.WeiCent.V||Para.WeiCent.R
            if isempty(dir([outdirCent,'VolWeiCent']))
                ClincFreq_CentPointForVolume_pre(INDIR,outdirCent,Para,spmjobdir);
            end
        end
    end
end
if Para.CentNum.R
    OutCentdir = [Outdir,filesep,'CenterNumber',filesep];
    mkdir(OutCentdir);
    disp('Now Calculate CenterNumber!')
    ClincFreq_CalCentNum(outdirCent,OutCentdir,SetUpPara,Para,AtlasTempU)
    % stat needed
    ClincFreq_CentNumStat(outdirCent,OutCentdir,SetUpPara,Para,AtlasTempU);
    disp('CenterNumber finished')
end
if Para.AffVolume.R
    OutAffVoldir = [Outdir,filesep,'AffectedVolume',filesep];
    mkdir(OutAffVoldir);
    disp('Now Calculate AffectedVolume!')
    ClincFreq_AffVolume(INDIR,OutAffVoldir,SetUpPara,Para,AtlasTempU)
    % stat needed
    ClincFreq_AffVolStat(INDIR,OutAffVoldir,SetUpPara,Para,AtlasTempU);
    disp('AffectedVolume finished')
end
if Para.WeiCent.V||Para.WeiCent.R
    OutWeiCentdir = [Outdir,filesep,'VolumeWeightedCenter',filesep];
    mkdir(OutWeiCentdir);
    disp('Now Calculate VolumeWeightedCenter!')
    ClincFreq_CalWeiCent(outdirCent,OutWeiCentdir,SetUpPara,Para,AtlasTempU,dbmask)
    % stat needed
    ClincFreq_WeiCentStat(outdirCent,OutWeiCentdir,SetUpPara,Para,AtlasTempU,dbmask);
    disp('VolumeWeightedCenter finished')
end

if Para.WeiLoc.Vol
    OutWeiLocdir = [Outdir,filesep,'VolumeWeightedLocation',filesep];
    mkdir(OutWeiLocdir);
    disp('Now Calculate Volume Weighted Location!')
    ClincFreq_CalClinWeiLocation(INDIR,OutWeiLocdir,SetUpPara,Para,1,[]);
    disp('Now Do the statistical analysis of Volume Weighted Location!')
    ClincFreq_ClinWeiLocationStat(INDIR,OutWeiLocdir,SetUpPara,Para,1,[]);
    disp('Volume Weighted Location finished!')
end
if Para.WeiLoc.Clin
    for icn = 1:ClincNum
        OutWeiLocClincdir = [Outdir,filesep,'ClincWeightedLocation_',ClincWei(icn).Name,filesep];
        mkdir(OutWeiLocClincdir)
        disp(['Now Calculate ',ClincWei(icn).Name,' Weighted Location!'])
        ClincFreq_CalClinWeiLocation(INDIR,OutWeiLocClincdir,SetUpPara,Para,0,ClincWei(icn));
        disp(['Now Do the statistical analysis of ',ClincWei(icn).Name,' Weighted Location!'])
        ClincFreq_ClinWeiLocationStat(INDIR,OutWeiLocClincdir,SetUpPara,Para,0,ClincWei(icn));
        disp([ClincWei(icn).Name,' Weighted Location finished!'])
    end
end

if Para.Other.Hist
    OutHistdir = [Outdir,filesep,'Histgram',filesep];
    mkdir(OutHistdir);
    ClincFreq_Hist(INDIR,INDIRBG,OutHistdir,voxelsize);
end
if Para.Other.D2V
    OutD2Vdir = [Outdir,filesep,'DistanceToVentricle',filesep];
    mkdir(OutD2Vdir);
    ClincFreq_D2V(INDIR,OutD2Vdir,voxelsize,[Outdir,filesep,'ReslicedAtlas',filesep,'AtlasU_Ventricle.nii'],[Outdir,filesep,'ReslicedAtlas',filesep,'AtlasU_Ventricle.mat']);
end

disp('Job finished!')

end