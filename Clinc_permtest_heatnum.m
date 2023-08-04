function [pval ClusNum] = Clinc_permtest_heatnum(DATtemp1,DATtemp2,NUM1,NUM2,Lab,V,VoxelDirPermcluster,Dmasked)
% Lab: 1 voxel, 2 ROI
pths = which('Clinc_permtest.m');
[pth,fil,nam] = fileparts(pths);
if Lab==1  %
    if isempty(dir([pth,filesep,'permtest_clust']))
        mkdir([pth,filesep,'permtest_clust']);
    else
        rmdir([pth,filesep,'permtest_clust'],'s');
        mkdir([pth,filesep,'permtest_clust']);
    end
end
% permclust_dir = [pth,filesep,'permtest_clust',filesep];
permclust_dir = [VoxelDirPermcluster,filesep];
% pernum = 5000;
pernum = 1000;
DATtempall = [sparse(DATtemp1),sparse(DATtemp2)];
TOTALNUM = NUM1+NUM2;
inds = find(sum(DATtemp1,2)+sum(DATtemp2,2));
DATu = DATtempall(inds,:);
origd = sum(DATtemp1(inds,:),2)-sum(DATtemp2(inds,:),2);
% deltV = zeros(length(inds),5000);
if ~isempty(inds)
    cutleng = 10000;
%     cutleng = 5000;
    cutnum = floor(length(inds)/cutleng);
    for i = 1:cutnum
        CUTpiece{i,1} = 1+cutleng*(i-1):cutleng*i;
    end
    if length(inds)-cutnum*cutleng>0
        cutnum = cutnum+1;
        CUTpiece{cutnum,1} = (cutnum-1)*cutleng+1:length(inds);
    end
    tic
    Pval = zeros(length(inds),1);
    for icut = 1:cutnum
        Pvalperm = zeros(length(CUTpiece{icut}),pernum);
        DATu2 = DATu(CUTpiece{icut},:);
        deltV = zeros(length(CUTpiece{icut}),pernum);
        parfor iperm = 1:pernum
            TOTALNUMp = randperm(TOTALNUM);
            G1 = TOTALNUMp(1:NUM1);
            G2 = TOTALNUMp(NUM1+1:TOTALNUM);
            deltV(:,iperm) = sum(DATu2(:,G1),2)-sum(DATu2(:,G2),2);
        end
        origDtemp = origd(CUTpiece{icut});
        [mu sig muci sigci] = normfit(deltV');
        Pval(CUTpiece{icut}) = normcdf(origDtemp',mu,sig);
        if Lab==1
            parfor iperm = 1:pernum
                origDtemp_P = deltV(:,iperm);
                Pvalperm(:,iperm) = normcdf(origDtemp_P',mu,sig);
            end
            save([permclust_dir,'Piece',sprintf('%05d',icut),'.mat'],'Pvalperm','CUTpiece','icut');
        end
        clear origDtemp deltV origDtemp_P Pvalperm
        toc
    end
    pval = zeros(size(DATtemp1,1),1);
    pval(inds) = Pval;
    if Lab==1
        tic
        %%
%         CutPermNum = 100;
%         CutPermR = floor(pernum/CutPermNum);
%         for icut = 1:CutPermR
%             CutPiece{icut,1} = [1:CutPermNum]+CutPermNum*(icut-1);
%         end
%         if CutPermNum*CutPermR<pernum
%             CutPermR = CutPermR+1;
%             CutPiece{CutPermR,1} = (CutPermNum*CutPermR+1):pernum;
%         end
%         for icut = 1:CutPermR
%             CutPtemp = CutPiece{icut};
%             for jcut = 1:length(CutPtemp)
%                 pvalperm = zeros(size(DATtemp1,1),1);
%                 Pvalp = getinfomat(cutnum,permclust_dir,inds,CutPtemp(jcut));
%                 pvalperm(inds) = Pvalp;
%                 pvaltemp(jcut,:) = pvalperm;
%                 pvalperm = pvalperm.*reshape(Dmasked,prod(V(1).dim),1);
%                 pvaltempIG(jcut,:) = pvalperm;
%             end
%             save([permclust_dir,'Permpiece',sprintf('%05d',icut),'.mat'],'pvaltemp','CutPtemp','pvaltempIG');
%             clear pvaltemp pvaltempIG;
%         end
%         parfor iperm = 1:pernum
%             
%         end
        %%
        try
            parfor iperm = 1:pernum
                pvalperm = zeros(size(DATtemp1,1),1);
                Pvalp = getinfomat(cutnum,permclust_dir,inds,iperm);
                
                pvalperm(inds) = Pvalp;
                PvalpermV = reshape(pvalperm,V(1).dim);
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05)+PvalpermV>(1-0.05);
                PvalpermV005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01)+PvalpermV>(1-0.01);
                PvalpermV001Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005)+PvalpermV>(1-0.005);
                PvalpermV0005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001)+PvalpermV>(1-0.001);
                PvalpermV0001Ct{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = PvalpermV>(1-0.05);
                PvalpermV005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.01);
                PvalpermV001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.005);
                PvalpermV0005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.001);
                PvalpermV0001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05);
                PvalpermV005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01);
                PvalpermV001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005);
                PvalpermV0005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001);
                PvalpermV0001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                
                %% Ignore Percent     
                PvalpermV = PvalpermV.*Dmasked;
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05)+PvalpermV>(1-0.05);
                IG_PvalpermV005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01)+PvalpermV>(1-0.01);
                IG_PvalpermV001Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005)+PvalpermV>(1-0.005);
                IG_PvalpermV0005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001)+PvalpermV>(1-0.001);
                IG_PvalpermV0001Ct{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = PvalpermV>(1-0.05);
                IG_PvalpermV005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.01);
                IG_PvalpermV001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.005);
                IG_PvalpermV0005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.001);
                IG_PvalpermV0001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05);
                IG_PvalpermV005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01);
                IG_PvalpermV001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005);
                IG_PvalpermV0005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001);
                IG_PvalpermV0001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                                
                %             toc
            end
            toc
        catch
            for iperm = 1:pernum
                pvalperm = zeros(size(DATtemp1,1),1);
                Pvalp = getinfomat(cutnum,permclust_dir,inds,iperm);
                
                pvalperm(inds) = Pvalp;
                PvalpermV = reshape(pvalperm,V(1).dim);
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05)+PvalpermV>(1-0.05);
                PvalpermV005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01)+PvalpermV>(1-0.01);
                PvalpermV001Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005)+PvalpermV>(1-0.005);
                PvalpermV0005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001)+PvalpermV>(1-0.001);
                PvalpermV0001Ct{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = PvalpermV>(1-0.05);
                PvalpermV005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.01);
                PvalpermV001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.005);
                PvalpermV0005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.001);
                PvalpermV0001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05);
                PvalpermV005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01);
                PvalpermV001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005);
                PvalpermV0005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001);
                PvalpermV0001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                
                %% Ignore Percent     
                PvalpermV = PvalpermV.*Dmasked;
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05)+PvalpermV>(1-0.05);
                IG_PvalpermV005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01)+PvalpermV>(1-0.01);
                IG_PvalpermV001Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005)+PvalpermV>(1-0.005);
                IG_PvalpermV0005Ct{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001)+PvalpermV>(1-0.001);
                IG_PvalpermV0001Ct{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = PvalpermV>(1-0.05);
                IG_PvalpermV005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.01);
                IG_PvalpermV001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.005);
                IG_PvalpermV0005Ctp{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = PvalpermV>(1-0.001);
                IG_PvalpermV0001Ctp{iperm,1} = getclusterinfo2(PvalpermVt);                
                
                PvalpermVt = (PvalpermV>0&PvalpermV<0.05);
                IG_PvalpermV005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.01);
                IG_PvalpermV001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.005);
                IG_PvalpermV0005Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
                PvalpermVt = (PvalpermV>0&PvalpermV<0.001);
                IG_PvalpermV0001Ctn{iperm,1} = getclusterinfo2(PvalpermVt);
            end
            toc
        end
        save([permclust_dir,'clusterinfo.mat'],'PvalpermV0001Ct','PvalpermV0005Ct','PvalpermV001Ct','PvalpermV005Ct',...
            'PvalpermV0001Ctp','PvalpermV0005Ctp','PvalpermV001Ctp','PvalpermV005Ctp',...
            'PvalpermV0001Ctn','PvalpermV0005Ctn','PvalpermV001Ctn','PvalpermV005Ctn',...
            'IG_PvalpermV0001Ct','IG_PvalpermV0005Ct','IG_PvalpermV001Ct','IG_PvalpermV005Ct',...
            'IG_PvalpermV0001Ctp','IG_PvalpermV0005Ctp','IG_PvalpermV001Ctp','IG_PvalpermV005Ctp',...
            'IG_PvalpermV0001Ctn','IG_PvalpermV0005Ctn','IG_PvalpermV001Ctn','IG_PvalpermV005Ctn')
        ClusNum.ClusterNumB0001 = PermClusterInfo(PvalpermV0001Ct);
        ClusNum.ClusterNumB0005 = PermClusterInfo(PvalpermV0005Ct);
        ClusNum.ClusterNumB001 = PermClusterInfo(PvalpermV001Ct);
        ClusNum.ClusterNumB005 = PermClusterInfo(PvalpermV005Ct);
        ClusNum.ClusterNumP0001 = PermClusterInfo(PvalpermV0001Ctp);
        ClusNum.ClusterNumP0005 = PermClusterInfo(PvalpermV0005Ctp);
        ClusNum.ClusterNumP001 = PermClusterInfo(PvalpermV001Ctp);
        ClusNum.ClusterNumP005 = PermClusterInfo(PvalpermV005Ctp);
        ClusNum.ClusterNumN0001 = PermClusterInfo(PvalpermV0001Ctn);
        ClusNum.ClusterNumN0005 = PermClusterInfo(PvalpermV0005Ctn);
        ClusNum.ClusterNumN001 = PermClusterInfo(PvalpermV001Ctn);
        ClusNum.ClusterNumN005 = PermClusterInfo(PvalpermV005Ctn);
        %% Ignore Percent
        ClusNum.IG_ClusterNumB0001 = PermClusterInfo(IG_PvalpermV0001Ct);
        ClusNum.IG_ClusterNumB0005 = PermClusterInfo(IG_PvalpermV0005Ct);
        ClusNum.IG_ClusterNumB001 = PermClusterInfo(IG_PvalpermV001Ct);
        ClusNum.IG_ClusterNumB005 = PermClusterInfo(IG_PvalpermV005Ct);
        ClusNum.IG_ClusterNumP0001 = PermClusterInfo(IG_PvalpermV0001Ctp);
        ClusNum.IG_ClusterNumP0005 = PermClusterInfo(IG_PvalpermV0005Ctp);
        ClusNum.IG_ClusterNumP001 = PermClusterInfo(IG_PvalpermV001Ctp);
        ClusNum.IG_ClusterNumP005 = PermClusterInfo(IG_PvalpermV005Ctp);
        ClusNum.IG_ClusterNumN0001 = PermClusterInfo(IG_PvalpermV0001Ctn);
        ClusNum.IG_ClusterNumN0005 = PermClusterInfo(IG_PvalpermV0005Ctn);
        ClusNum.IG_ClusterNumN001 = PermClusterInfo(IG_PvalpermV001Ctn);
        ClusNum.IG_ClusterNumN005 = PermClusterInfo(IG_PvalpermV005Ctn);
    else
        ClusNum = [];
    end
    % for i = 1:length(INDs)
    %     [mu sig muci sigci] = normfit(deltG1G2PermTotal(i,:));
    %     Pval(i) = normcdf(deltG1G2(INDs(i)),mu,sig);
    % end
else
    pval = zeros(size(DATtemp1,1),1);
    ClusNum = [];
end
end