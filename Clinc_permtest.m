function [pval ClusNum] = Clinc_permtest(DATtemp1,DATtemp2,NUM1,NUM2,Lab,V,VoxelDirPermcluster)
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
        try
            parfor iperm = 1:pernum
                pvalperm = zeros(size(DATtemp1,1),1);
                Pvalp = getinfomat(cutnum,permclust_dir,inds,iperm);
                
                pvalperm(inds) = Pvalp;
                PvalpermV = reshape(pvalperm,V(1).dim);
                PvalpermV005 = (PvalpermV>0&PvalpermV<0.05)+PvalpermV>(1-0.05);
                PvalpermV001 = (PvalpermV>0&PvalpermV<0.01)+PvalpermV>(1-0.01);
                PvalpermV0005 = (PvalpermV>0&PvalpermV<0.005)+PvalpermV>(1-0.005);
                PvalpermV0001 = (PvalpermV>0&PvalpermV<0.001)+PvalpermV>(1-0.001);
                
                PvalpermV005p = PvalpermV>(1-0.05);
                PvalpermV001p = PvalpermV>(1-0.01);
                PvalpermV0005p = PvalpermV>(1-0.005);
                PvalpermV0001p = PvalpermV>(1-0.001);
                
                PvalpermV005n = (PvalpermV>0&PvalpermV<0.05);
                PvalpermV001n = (PvalpermV>0&PvalpermV<0.01);
                PvalpermV0005n = (PvalpermV>0&PvalpermV<0.005);
                PvalpermV0001n = (PvalpermV>0&PvalpermV<0.001);
                
                [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
                    getclusterinfo(PvalpermV005,PvalpermV001,PvalpermV0005,PvalpermV0001);
                PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
                PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
                PvalpermV001Ct{iperm,1} = PvalpermV001C;
                PvalpermV005Ct{iperm,1} = PvalpermV005C;
                
                [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
                    getclusterinfo(PvalpermV005p,PvalpermV001p,PvalpermV0005p,PvalpermV0001p);
                PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
                PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
                PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
                PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
                
                [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
                    getclusterinfo(PvalpermV005n,PvalpermV001n,PvalpermV0005n,PvalpermV0001n);
                PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
                PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
                PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
                PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;
                %             toc
            end
            toc
        catch
            for iperm = 1:pernum
                pvalperm = zeros(size(DATtemp1,1),1);
                Pvalp = getinfomat(cutnum,permclust_dir,inds,iperm);
                
                pvalperm(inds) = Pvalp;
                PvalpermV = reshape(pvalperm,V(1).dim);
                PvalpermV005 = (PvalpermV>0&PvalpermV<0.05)+PvalpermV>(1-0.05);
                PvalpermV001 = (PvalpermV>0&PvalpermV<0.01)+PvalpermV>(1-0.01);
                PvalpermV0005 = (PvalpermV>0&PvalpermV<0.005)+PvalpermV>(1-0.005);
                PvalpermV0001 = (PvalpermV>0&PvalpermV<0.001)+PvalpermV>(1-0.001);
                
                PvalpermV005p = PvalpermV>(1-0.05);
                PvalpermV001p = PvalpermV>(1-0.01);
                PvalpermV0005p = PvalpermV>(1-0.005);
                PvalpermV0001p = PvalpermV>(1-0.001);
                
                PvalpermV005n = (PvalpermV>0&PvalpermV<0.05);
                PvalpermV001n = (PvalpermV>0&PvalpermV<0.01);
                PvalpermV0005n = (PvalpermV>0&PvalpermV<0.005);
                PvalpermV0001n = (PvalpermV>0&PvalpermV<0.001);
                
                [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
                    getclusterinfo(PvalpermV005,PvalpermV001,PvalpermV0005,PvalpermV0001);
                PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
                PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
                PvalpermV001Ct{iperm,1} = PvalpermV001C;
                PvalpermV005Ct{iperm,1} = PvalpermV005C;
                
                [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
                    getclusterinfo(PvalpermV005p,PvalpermV001p,PvalpermV0005p,PvalpermV0001p);
                PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
                PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
                PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
                PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
                
                [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
                    getclusterinfo(PvalpermV005n,PvalpermV001n,PvalpermV0005n,PvalpermV0001n);
                PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
                PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
                PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
                PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;
                %             toc
            end
            toc
        end
        save([permclust_dir,'clusterinfo.mat'],'PvalpermV0001Ct','PvalpermV0005Ct','PvalpermV001Ct','PvalpermV005Ct',...
            'PvalpermV0001Ctp','PvalpermV0005Ctp','PvalpermV001Ctp','PvalpermV005Ctp',...
            'PvalpermV0001Ctn','PvalpermV0005Ctn','PvalpermV001Ctn','PvalpermV005Ctn')
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