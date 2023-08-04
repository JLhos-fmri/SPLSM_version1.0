function ExtClustInfo_CWL(CWLclusterDir,C,j,Dmasked,OWLDir)
% save([OWLclus{i},filesep,'Piece',sprintf('%05d',iseg),'.mat'],...
%     'Pval_media_perm','Pval_mean_perm','Pval_std_perm','segind3','iseg','DATIND','V','Dmasked');
piecedir = [CWLclusterDir,filesep,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),filesep];
pieces_files = dir([CWLclusterDir,filesep,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),filesep,'Piece*.mat']);
load([CWLclusterDir,filesep,'permnums.mat']);
for iperm = 1:permnums
    for i = 1:length(pieces_files)
        load([piecedir,filesep,pieces_files(i).name]);
        PvalClust_media(segind3{iseg}) = Pval_media_perm(iperm,:);
        PvalClust_mean(segind3{iseg}) = Pval_mean_perm(iperm,:);
        PvalClust_std(segind3{iseg}) = Pval_std_perm(iperm,:);
    end
    PvalC_media = zeros(V(1).dim);
    PvalC_media(DATIND) = PvalClust_media;
    PvalC_media_Pboth05 = (PvalC_media>0.95)+(PvalC_media<0.05&PvalC_media>0);
    PvalC_media_Pboth01 = (PvalC_media>0.99)+(PvalC_media<0.01&PvalC_media>0);
    PvalC_media_Pboth005 = (PvalC_media>0.995)+(PvalC_media<0.005&PvalC_media>0);
    PvalC_media_Pboth001 = (PvalC_media>0.999)+(PvalC_media<0.001&PvalC_media>0);
    
    PvalC_media_Ppos05 = (PvalC_media>0.95);
    PvalC_media_Ppos01 = (PvalC_media>0.99);
    PvalC_media_Ppos005 = (PvalC_media>0.995);
    PvalC_media_Ppos001 = (PvalC_media>0.999);
    
    PvalC_media_Pneg05 = (PvalC_media<0.05&PvalC_media>0);
    PvalC_media_Pneg01 = (PvalC_media<0.01&PvalC_media>0);
    PvalC_media_Pneg005 = (PvalC_media<0.005&PvalC_media>0);
    PvalC_media_Pneg001 = (PvalC_media<0.001&PvalC_media>0);    
    
    [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
        getclusterinfo(PvalC_media_Pboth05,PvalC_media_Pboth01,PvalC_media_Pboth005,PvalC_media_Pboth001);
    media_PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
    media_PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
    media_PvalpermV001Ct{iperm,1} = PvalpermV001C;
    media_PvalpermV005Ct{iperm,1} = PvalpermV005C;    
    
    [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
        getclusterinfo(PvalC_media_Ppos05,PvalC_media_Ppos01,PvalC_media_Ppos005,PvalC_media_Ppos001);
    media_PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
    media_PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
    media_PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
    media_PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
    
    [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
        getclusterinfo(PvalC_media_Pneg05,PvalC_media_Pneg01,PvalC_media_Pneg005,PvalC_media_Pneg001);
    media_PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
    media_PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
    media_PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
    media_PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;    
    
    %%    
    PvalC_mean = zeros(V(1).dim);
    PvalC_mean(DATIND) = PvalClust_mean;
    PvalC_mean_Pboth05 = (PvalC_mean>0.95)+(PvalC_mean<0.05&PvalC_mean>0);
    PvalC_mean_Pboth01 = (PvalC_mean>0.99)+(PvalC_mean<0.01&PvalC_mean>0);
    PvalC_mean_Pboth005 = (PvalC_mean>0.995)+(PvalC_mean<0.005&PvalC_mean>0);
    PvalC_mean_Pboth001 = (PvalC_mean>0.999)+(PvalC_mean<0.001&PvalC_mean>0);
    
    PvalC_mean_Ppos05 = (PvalC_mean>0.95);
    PvalC_mean_Ppos01 = (PvalC_mean>0.99);
    PvalC_mean_Ppos005 = (PvalC_mean>0.995);
    PvalC_mean_Ppos001 = (PvalC_mean>0.999);
    
    PvalC_mean_Pneg05 = (PvalC_mean<0.05&PvalC_mean>0);
    PvalC_mean_Pneg01 = (PvalC_mean<0.01&PvalC_mean>0);
    PvalC_mean_Pneg005 = (PvalC_mean<0.005&PvalC_mean>0);
    PvalC_mean_Pneg001 = (PvalC_mean<0.001&PvalC_mean>0);    
    
    [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
        getclusterinfo(PvalC_mean_Pboth05,PvalC_mean_Pboth01,PvalC_mean_Pboth005,PvalC_mean_Pboth001);
    mean_PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
    mean_PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
    mean_PvalpermV001Ct{iperm,1} = PvalpermV001C;
    mean_PvalpermV005Ct{iperm,1} = PvalpermV005C;    
    
    [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
        getclusterinfo(PvalC_mean_Ppos05,PvalC_mean_Ppos01,PvalC_mean_Ppos005,PvalC_mean_Ppos001);
    mean_PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
    mean_PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
    mean_PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
    mean_PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
    
    [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
        getclusterinfo(PvalC_mean_Pneg05,PvalC_mean_Pneg01,PvalC_mean_Pneg005,PvalC_mean_Pneg001);
    mean_PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
    mean_PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
    mean_PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
    mean_PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;
    
    %%    
    PvalC_std = zeros(V(1).dim);
    PvalC_std(DATIND) = PvalClust_std;
    PvalC_std_Pboth05 = (PvalC_std>0.95)+(PvalC_std<0.05&PvalC_std>0);
    PvalC_std_Pboth01 = (PvalC_std>0.99)+(PvalC_std<0.01&PvalC_std>0);
    PvalC_std_Pboth005 = (PvalC_std>0.995)+(PvalC_std<0.005&PvalC_std>0);
    PvalC_std_Pboth001 = (PvalC_std>0.999)+(PvalC_std<0.001&PvalC_std>0);
    
    PvalC_std_Ppos05 = (PvalC_std>0.95);
    PvalC_std_Ppos01 = (PvalC_std>0.99);
    PvalC_std_Ppos005 = (PvalC_std>0.995);
    PvalC_std_Ppos001 = (PvalC_std>0.999);
    
    PvalC_std_Pneg05 = (PvalC_std<0.05&PvalC_std>0);
    PvalC_std_Pneg01 = (PvalC_std<0.01&PvalC_std>0);
    PvalC_std_Pneg005 = (PvalC_std<0.005&PvalC_std>0);
    PvalC_std_Pneg001 = (PvalC_std<0.001&PvalC_std>0);    
    
    [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
        getclusterinfo(PvalC_std_Pboth05,PvalC_std_Pboth01,PvalC_std_Pboth005,PvalC_std_Pboth001);
    std_PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
    std_PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
    std_PvalpermV001Ct{iperm,1} = PvalpermV001C;
    std_PvalpermV005Ct{iperm,1} = PvalpermV005C;    
    
    [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
        getclusterinfo(PvalC_std_Ppos05,PvalC_std_Ppos01,PvalC_std_Ppos005,PvalC_std_Ppos001);
    std_PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
    std_PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
    std_PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
    std_PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
    
    [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
        getclusterinfo(PvalC_std_Pneg05,PvalC_std_Pneg01,PvalC_std_Pneg005,PvalC_std_Pneg001);
    std_PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
    std_PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
    std_PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
    std_PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;
    
    %% IG: 
    
    PvalC_media = zeros(V(1).dim);
    PvalC_media(DATIND) = PvalClust_media;
    PvalC_media_Pboth05 = (PvalC_media>0.95)+(PvalC_media<0.05&PvalC_media>0);
    PvalC_media_Pboth01 = (PvalC_media>0.99)+(PvalC_media<0.01&PvalC_media>0);
    PvalC_media_Pboth005 = (PvalC_media>0.995)+(PvalC_media<0.005&PvalC_media>0);
    PvalC_media_Pboth001 = (PvalC_media>0.999)+(PvalC_media<0.001&PvalC_media>0);
    
    PvalC_media_Ppos05 = (PvalC_media>0.95);
    PvalC_media_Ppos01 = (PvalC_media>0.99);
    PvalC_media_Ppos005 = (PvalC_media>0.995);
    PvalC_media_Ppos001 = (PvalC_media>0.999);
    
    PvalC_media_Pneg05 = (PvalC_media<0.05&PvalC_media>0);
    PvalC_media_Pneg01 = (PvalC_media<0.01&PvalC_media>0);
    PvalC_media_Pneg005 = (PvalC_media<0.005&PvalC_media>0);
    PvalC_media_Pneg001 = (PvalC_media<0.001&PvalC_media>0);    
    
    [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
        getclusterinfo(PvalC_media_Pboth05.*Dmasked,PvalC_media_Pboth01.*Dmasked,PvalC_media_Pboth005.*Dmasked,PvalC_media_Pboth001.*Dmasked);
    IG_media_PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
    IG_media_PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
    IG_media_PvalpermV001Ct{iperm,1} = PvalpermV001C;
    IG_media_PvalpermV005Ct{iperm,1} = PvalpermV005C;    
    
    [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
        getclusterinfo(PvalC_media_Ppos05.*Dmasked,PvalC_media_Ppos01.*Dmasked,PvalC_media_Ppos005.*Dmasked,PvalC_media_Ppos001.*Dmasked);
    IG_media_PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
    IG_media_PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
    IG_media_PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
    IG_media_PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
    
    [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
        getclusterinfo(PvalC_media_Pneg05.*Dmasked,PvalC_media_Pneg01.*Dmasked,PvalC_media_Pneg005.*Dmasked,PvalC_media_Pneg001.*Dmasked);
    IG_media_PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
    IG_media_PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
    IG_media_PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
    IG_media_PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;    
    
    %%    
    PvalC_mean = zeros(V(1).dim);
    PvalC_mean(DATIND) = PvalClust_mean;
    PvalC_mean_Pboth05 = (PvalC_mean>0.95)+(PvalC_mean<0.05&PvalC_mean>0);
    PvalC_mean_Pboth01 = (PvalC_mean>0.99)+(PvalC_mean<0.01&PvalC_mean>0);
    PvalC_mean_Pboth005 = (PvalC_mean>0.995)+(PvalC_mean<0.005&PvalC_mean>0);
    PvalC_mean_Pboth001 = (PvalC_mean>0.999)+(PvalC_mean<0.001&PvalC_mean>0);
    
    PvalC_mean_Ppos05 = (PvalC_mean>0.95);
    PvalC_mean_Ppos01 = (PvalC_mean>0.99);
    PvalC_mean_Ppos005 = (PvalC_mean>0.995);
    PvalC_mean_Ppos001 = (PvalC_mean>0.999);
    
    PvalC_mean_Pneg05 = (PvalC_mean<0.05&PvalC_mean>0);
    PvalC_mean_Pneg01 = (PvalC_mean<0.01&PvalC_mean>0);
    PvalC_mean_Pneg005 = (PvalC_mean<0.005&PvalC_mean>0);
    PvalC_mean_Pneg001 = (PvalC_mean<0.001&PvalC_mean>0);    
    
    [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
        getclusterinfo(PvalC_mean_Pboth05.*Dmasked,PvalC_mean_Pboth01.*Dmasked,PvalC_mean_Pboth005.*Dmasked,PvalC_mean_Pboth001.*Dmasked);
    IG_mean_PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
    IG_mean_PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
    IG_mean_PvalpermV001Ct{iperm,1} = PvalpermV001C;
    IG_mean_PvalpermV005Ct{iperm,1} = PvalpermV005C;    
    
    [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
        getclusterinfo(PvalC_mean_Ppos05.*Dmasked,PvalC_mean_Ppos01.*Dmasked,PvalC_mean_Ppos005.*Dmasked,PvalC_mean_Ppos001.*Dmasked);
    IG_mean_PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
    IG_mean_PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
    IG_mean_PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
    IG_mean_PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
    
    [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
        getclusterinfo(PvalC_mean_Pneg05.*Dmasked,PvalC_mean_Pneg01.*Dmasked,PvalC_mean_Pneg005.*Dmasked,PvalC_mean_Pneg001.*Dmasked);
    IG_mean_PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
    IG_mean_PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
    IG_mean_PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
    IG_mean_PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;
    
    %%    
    PvalC_std = zeros(V(1).dim);
    PvalC_std(DATIND) = PvalClust_std;
    PvalC_std_Pboth05 = (PvalC_std>0.95)+(PvalC_std<0.05&PvalC_std>0);
    PvalC_std_Pboth01 = (PvalC_std>0.99)+(PvalC_std<0.01&PvalC_std>0);
    PvalC_std_Pboth005 = (PvalC_std>0.995)+(PvalC_std<0.005&PvalC_std>0);
    PvalC_std_Pboth001 = (PvalC_std>0.999)+(PvalC_std<0.001&PvalC_std>0);
    
    PvalC_std_Ppos05 = (PvalC_std>0.95);
    PvalC_std_Ppos01 = (PvalC_std>0.99);
    PvalC_std_Ppos005 = (PvalC_std>0.995);
    PvalC_std_Ppos001 = (PvalC_std>0.999);
    
    PvalC_std_Pneg05 = (PvalC_std<0.05&PvalC_std>0);
    PvalC_std_Pneg01 = (PvalC_std<0.01&PvalC_std>0);
    PvalC_std_Pneg005 = (PvalC_std<0.005&PvalC_std>0);
    PvalC_std_Pneg001 = (PvalC_std<0.001&PvalC_std>0);    
    
    [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = ...
        getclusterinfo(PvalC_std_Pboth05.*Dmasked,PvalC_std_Pboth01.*Dmasked,PvalC_std_Pboth005.*Dmasked,PvalC_std_Pboth001.*Dmasked);
    IG_std_PvalpermV0001Ct{iperm,1} = PvalpermV0001C;
    IG_std_PvalpermV0005Ct{iperm,1} = PvalpermV0005C;
    IG_std_PvalpermV001Ct{iperm,1} = PvalpermV001C;
    IG_std_PvalpermV005Ct{iperm,1} = PvalpermV005C;    
    
    [PvalpermV0001Cp,PvalpermV0005Cp,PvalpermV001Cp,PvalpermV005Cp] = ...
        getclusterinfo(PvalC_std_Ppos05.*Dmasked,PvalC_std_Ppos01.*Dmasked,PvalC_std_Ppos005.*Dmasked,PvalC_std_Ppos001.*Dmasked);
    IG_std_PvalpermV0001Ctp{iperm,1} = PvalpermV0001Cp;
    IG_std_PvalpermV0005Ctp{iperm,1} = PvalpermV0005Cp;
    IG_std_PvalpermV001Ctp{iperm,1} = PvalpermV001Cp;
    IG_std_PvalpermV005Ctp{iperm,1} = PvalpermV005Cp;
    
    [PvalpermV0001Cn,PvalpermV0005Cn,PvalpermV001Cn,PvalpermV005Cn] = ...
        getclusterinfo(PvalC_std_Pneg05.*Dmasked,PvalC_std_Pneg01.*Dmasked,PvalC_std_Pneg005.*Dmasked,PvalC_std_Pneg001.*Dmasked);
    IG_std_PvalpermV0001Ctn{iperm,1} = PvalpermV0001Cn;
    IG_std_PvalpermV0005Ctn{iperm,1} = PvalpermV0005Cn;
    IG_std_PvalpermV001Ctn{iperm,1} = PvalpermV001Cn;
    IG_std_PvalpermV005Ctn{iperm,1} = PvalpermV005Cn;
    
end
save([CWLclusterDir,filesep,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),filesep,'clusterinfo.mat'],...
    'media_PvalpermV0001Ct','media_PvalpermV0005Ct','media_PvalpermV005Ct','media_PvalpermV001Ct',...
    'media_PvalpermV0001Ctp','media_PvalpermV0005Ctp','media_PvalpermV005Ctp','media_PvalpermV001Ctp',...
    'media_PvalpermV0001Ctn','media_PvalpermV0005Ctn','media_PvalpermV005Ctn','media_PvalpermV001Ctn',...
    'mean_PvalpermV0001Ct','mean_PvalpermV0005Ct','mean_PvalpermV005Ct','mean_PvalpermV001Ct',...
    'mean_PvalpermV0001Ctp','mean_PvalpermV0005Ctp','mean_PvalpermV005Ctp','mean_PvalpermV001Ctp',...
    'mean_PvalpermV0001Ctn','mean_PvalpermV0005Ctn','mean_PvalpermV005Ctn','mean_PvalpermV001Ctn',...
    'std_PvalpermV0001Ct','std_PvalpermV0005Ct','std_PvalpermV005Ct','std_PvalpermV001Ct',...
    'std_PvalpermV0001Ctp','std_PvalpermV0005Ctp','std_PvalpermV005Ctp','std_PvalpermV001Ctp',...
    'std_PvalpermV0001Ctn','std_PvalpermV0005Ctn','std_PvalpermV005Ctn','std_PvalpermV001Ctn',...
    'IG_media_PvalpermV0001Ct','IG_media_PvalpermV0005Ct','IG_media_PvalpermV005Ct','IG_media_PvalpermV001Ct',...
    'IG_media_PvalpermV0001Ctp','IG_media_PvalpermV0005Ctp','IG_media_PvalpermV005Ctp','IG_media_PvalpermV001Ctp',...
    'IG_media_PvalpermV0001Ctn','IG_media_PvalpermV0005Ctn','IG_media_PvalpermV005Ctn','IG_media_PvalpermV001Ctn',...
    'IG_mean_PvalpermV0001Ct','IG_mean_PvalpermV0005Ct','IG_mean_PvalpermV005Ct','IG_mean_PvalpermV001Ct',...
    'IG_mean_PvalpermV0001Ctp','IG_mean_PvalpermV0005Ctp','IG_mean_PvalpermV005Ctp','IG_mean_PvalpermV001Ctp',...
    'IG_mean_PvalpermV0001Ctn','IG_mean_PvalpermV0005Ctn','IG_mean_PvalpermV005Ctn','IG_mean_PvalpermV001Ctn',...
    'IG_std_PvalpermV0001Ct','IG_std_PvalpermV0005Ct','IG_std_PvalpermV005Ct','IG_std_PvalpermV001Ct',...
    'IG_std_PvalpermV0001Ctp','IG_std_PvalpermV0005Ctp','IG_std_PvalpermV005Ctp','IG_std_PvalpermV001Ctp',...
    'IG_std_PvalpermV0001Ctn','IG_std_PvalpermV0005Ctn','IG_std_PvalpermV005Ctn','IG_std_PvalpermV001Ctn');

ClusNum.media_ClusterNumB0001 = PermClusterInfo(media_PvalpermV0001Ct);
ClusNum.media_ClusterNumB0005 = PermClusterInfo(media_PvalpermV0005Ct);
ClusNum.media_ClusterNumB001 = PermClusterInfo(media_PvalpermV001Ct);
ClusNum.media_ClusterNumB005 = PermClusterInfo(media_PvalpermV005Ct);
ClusNum.media_ClusterNumP0001 = PermClusterInfo(media_PvalpermV0001Ctp);
ClusNum.media_ClusterNumP0005 = PermClusterInfo(media_PvalpermV0005Ctp);
ClusNum.media_ClusterNumP001 = PermClusterInfo(media_PvalpermV001Ctp);
ClusNum.media_ClusterNumP005 = PermClusterInfo(media_PvalpermV005Ctp);
ClusNum.media_ClusterNumN0001 = PermClusterInfo(media_PvalpermV0001Ctn);
ClusNum.media_ClusterNumN0005 = PermClusterInfo(media_PvalpermV0005Ctn);
ClusNum.media_ClusterNumN001 = PermClusterInfo(media_PvalpermV001Ctn);
ClusNum.media_ClusterNumN005 = PermClusterInfo(media_PvalpermV005Ctn);
%% Ignore Percent
ClusNum.IG_media_ClusterNumB0001 = PermClusterInfo(IG_media_PvalpermV0001Ct);
ClusNum.IG_media_ClusterNumB0005 = PermClusterInfo(IG_media_PvalpermV0005Ct);
ClusNum.IG_media_ClusterNumB001 = PermClusterInfo(IG_media_PvalpermV001Ct);
ClusNum.IG_media_ClusterNumB005 = PermClusterInfo(IG_media_PvalpermV005Ct);
ClusNum.IG_media_ClusterNumP0001 = PermClusterInfo(IG_media_PvalpermV0001Ctp);
ClusNum.IG_media_ClusterNumP0005 = PermClusterInfo(IG_media_PvalpermV0005Ctp);
ClusNum.IG_media_ClusterNumP001 = PermClusterInfo(IG_media_PvalpermV001Ctp);
ClusNum.IG_media_ClusterNumP005 = PermClusterInfo(IG_media_PvalpermV005Ctp);
ClusNum.IG_media_ClusterNumN0001 = PermClusterInfo(IG_media_PvalpermV0001Ctn);
ClusNum.IG_media_ClusterNumN0005 = PermClusterInfo(IG_media_PvalpermV0005Ctn);
ClusNum.IG_media_ClusterNumN001 = PermClusterInfo(IG_media_PvalpermV001Ctn);
ClusNum.IG_media_ClusterNumN005 = PermClusterInfo(IG_media_PvalpermV005Ctn);
%%
%%
ClusNum.mean_ClusterNumB0001 = PermClusterInfo(mean_PvalpermV0001Ct);
ClusNum.mean_ClusterNumB0005 = PermClusterInfo(mean_PvalpermV0005Ct);
ClusNum.mean_ClusterNumB001 = PermClusterInfo(mean_PvalpermV001Ct);
ClusNum.mean_ClusterNumB005 = PermClusterInfo(mean_PvalpermV005Ct);
ClusNum.mean_ClusterNumP0001 = PermClusterInfo(mean_PvalpermV0001Ctp);
ClusNum.mean_ClusterNumP0005 = PermClusterInfo(mean_PvalpermV0005Ctp);
ClusNum.mean_ClusterNumP001 = PermClusterInfo(mean_PvalpermV001Ctp);
ClusNum.mean_ClusterNumP005 = PermClusterInfo(mean_PvalpermV005Ctp);
ClusNum.mean_ClusterNumN0001 = PermClusterInfo(mean_PvalpermV0001Ctn);
ClusNum.mean_ClusterNumN0005 = PermClusterInfo(mean_PvalpermV0005Ctn);
ClusNum.mean_ClusterNumN001 = PermClusterInfo(mean_PvalpermV001Ctn);
ClusNum.mean_ClusterNumN005 = PermClusterInfo(mean_PvalpermV005Ctn);
%% Ignore Percent
ClusNum.IG_mean_ClusterNumB0001 = PermClusterInfo(IG_mean_PvalpermV0001Ct);
ClusNum.IG_mean_ClusterNumB0005 = PermClusterInfo(IG_mean_PvalpermV0005Ct);
ClusNum.IG_mean_ClusterNumB001 = PermClusterInfo(IG_mean_PvalpermV001Ct);
ClusNum.IG_mean_ClusterNumB005 = PermClusterInfo(IG_mean_PvalpermV005Ct);
ClusNum.IG_mean_ClusterNumP0001 = PermClusterInfo(IG_mean_PvalpermV0001Ctp);
ClusNum.IG_mean_ClusterNumP0005 = PermClusterInfo(IG_mean_PvalpermV0005Ctp);
ClusNum.IG_mean_ClusterNumP001 = PermClusterInfo(IG_mean_PvalpermV001Ctp);
ClusNum.IG_mean_ClusterNumP005 = PermClusterInfo(IG_mean_PvalpermV005Ctp);
ClusNum.IG_mean_ClusterNumN0001 = PermClusterInfo(IG_mean_PvalpermV0001Ctn);
ClusNum.IG_mean_ClusterNumN0005 = PermClusterInfo(IG_mean_PvalpermV0005Ctn);
ClusNum.IG_mean_ClusterNumN001 = PermClusterInfo(IG_mean_PvalpermV001Ctn);
ClusNum.IG_mean_ClusterNumN005 = PermClusterInfo(IG_mean_PvalpermV005Ctn);
%%
%%

ClusNum.std_ClusterNumB0001 = PermClusterInfo(std_PvalpermV0001Ct);
ClusNum.std_ClusterNumB0005 = PermClusterInfo(std_PvalpermV0005Ct);
ClusNum.std_ClusterNumB001 = PermClusterInfo(std_PvalpermV001Ct);
ClusNum.std_ClusterNumB005 = PermClusterInfo(std_PvalpermV005Ct);
ClusNum.std_ClusterNumP0001 = PermClusterInfo(std_PvalpermV0001Ctp);
ClusNum.std_ClusterNumP0005 = PermClusterInfo(std_PvalpermV0005Ctp);
ClusNum.std_ClusterNumP001 = PermClusterInfo(std_PvalpermV001Ctp);
ClusNum.std_ClusterNumP005 = PermClusterInfo(std_PvalpermV005Ctp);
ClusNum.std_ClusterNumN0001 = PermClusterInfo(std_PvalpermV0001Ctn);
ClusNum.std_ClusterNumN0005 = PermClusterInfo(std_PvalpermV0005Ctn);
ClusNum.std_ClusterNumN001 = PermClusterInfo(std_PvalpermV001Ctn);
ClusNum.std_ClusterNumN005 = PermClusterInfo(std_PvalpermV005Ctn);
%% Ignore Percent
ClusNum.IG_std_ClusterNumB0001 = PermClusterInfo(IG_std_PvalpermV0001Ct);
ClusNum.IG_std_ClusterNumB0005 = PermClusterInfo(IG_std_PvalpermV0005Ct);
ClusNum.IG_std_ClusterNumB001 = PermClusterInfo(IG_std_PvalpermV001Ct);
ClusNum.IG_std_ClusterNumB005 = PermClusterInfo(IG_std_PvalpermV005Ct);
ClusNum.IG_std_ClusterNumP0001 = PermClusterInfo(IG_std_PvalpermV0001Ctp);
ClusNum.IG_std_ClusterNumP0005 = PermClusterInfo(IG_std_PvalpermV0005Ctp);
ClusNum.IG_std_ClusterNumP001 = PermClusterInfo(IG_std_PvalpermV001Ctp);
ClusNum.IG_std_ClusterNumP005 = PermClusterInfo(IG_std_PvalpermV005Ctp);
ClusNum.IG_std_ClusterNumN0001 = PermClusterInfo(IG_std_PvalpermV0001Ctn);
ClusNum.IG_std_ClusterNumN0005 = PermClusterInfo(IG_std_PvalpermV0005Ctn);
ClusNum.IG_std_ClusterNumN001 = PermClusterInfo(IG_std_PvalpermV001Ctn);
ClusNum.IG_std_ClusterNumN005 = PermClusterInfo(IG_std_PvalpermV005Ctn);

% piecedir = [CWLclusterDir,filesep,'Group_',sprintf('%03d',C(j,1)),'vs_Group',sprintf('%03d',C(j,2)),filesep];
save([OWLDir,filesep,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'clusnum.mat'],'ClusNum')

end