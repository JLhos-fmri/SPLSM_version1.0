function ClusNum = VLSM_clusterPermutation(voxdir)
A = load([voxdir,filesep,'ForPerm.mat']);
permnum = 1000;
Dind = full(A.Dind);
GROUP1 = A.GROUP1;
GROUP2 = A.GROUP2;
COVterm = A.COVterm;
VALt = A.VALt;
indexs = A.indexs;
V = A.V;
tic
for iperm = 1:permnum
    ord = round(rand(1,size(Dind,2))*size(Dind,2));    
    ord(ord==0) = 1;
    ord(ord>size(Dind,2)) = size(Dind,2);
    t = zeros(length(indexs),1);
    p = ones(length(indexs),1)*0.5;
    Z = zeros(length(indexs),1);
    P = ones(length(indexs),1)*0.5;
    parfor j = 1:length(indexs)
        DAT = Dind(j,ord);
        [t(j),p(j),Z(j),P(j)] = VLSMsubfunc(DAT,GROUP1,GROUP2,COVterm(ord,:),VALt(ord,:));
    end
    
    Tmap = zeros(V(1).dim);
    Pmap = zeros(V(1).dim);
    Zmap = zeros(V(1).dim);
    Tmap(indexs) = t;
    Pmap(indexs) = p;
    Zmap(indexs) = Z;
    
    
%     pvalperm(inds) = Pvalp;
    PvalpermV = Pmap;
    PvalpermV005 = (PvalpermV<0.05&Tmap~=0);
    PvalpermV001 = (PvalpermV<0.01&Tmap~=0);
    PvalpermV0005 = (PvalpermV<0.005&Tmap~=0);
    PvalpermV0001 = (PvalpermV<0.001&Tmap~=0);
    
    PvalpermV005p = (PvalpermV<0.05&Tmap>0);
    PvalpermV001p = (PvalpermV<0.01&Tmap>0);
    PvalpermV0005p = (PvalpermV<0.005&Tmap>0);
    PvalpermV0001p = (PvalpermV<0.001&Tmap>0);
    
    PvalpermV005n = (PvalpermV<0.05&Tmap<0);
    PvalpermV001n = (PvalpermV<0.01&Tmap<0);
    PvalpermV0005n = (PvalpermV<0.005&Tmap<0);
    PvalpermV0001n = (PvalpermV<0.001&Tmap<0);
    
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
    toc
end
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

end