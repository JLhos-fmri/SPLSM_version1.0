function ClincFreq_Hist(INDIR,INDIRT1,OutHistdir,voxelsize)
files2 = dir([INDIR,'*.nii']);
files = dir([INDIRT1,'*.nii']);
xlsout{1,1} = 'pat';
xlsout{1,2} = 'voxel number';
xlsout{1,3} = 'volume(ml)';
xlsout{1,4} = 'min';
xlsout{1,5} = 'max';
xlsout{1,6} = 'mean';
xlsout{1,7} = 'std';
xlsout{1,8} = 'meadian';
xlsout{1,9} = 'Q1(25%)';
xlsout{1,10} = 'Q3(75%)';
xlsout{1,11} = 'kustosis';
xlsout{1,12} = 'skewness';
for i = 1:length(files)
    [v1 dat1] = Dynamic_read_dir_NIFTI([INDIRT1,files(i).name]);
    [v2 dat2] = Dynamic_read_dir_NIFTI([INDIR,files2(i).name]);

    voxnum = nnz(dat2);
    vsize = [1,1,1]*voxelsize;
    vols = prod(vsize)*voxnum/1000;
    indexs = find(dat2);
    datinroi = dat1(indexs);
    datinroi(isnan(datinroi)) = [];
    datinroi(isinf(datinroi)) = [];
    
    minv = min(datinroi);
    maxv = max(datinroi);
    meanv = mean(datinroi);
    stdv = mean(datinroi);
    medianv = median(datinroi);
    datsort = sort(datinroi);
    Q1v = datsort(round(length(datsort)/4));
    Q3v = datsort(round(length(datsort)/4*3));
    KUR = kurtosis(datinroi);
    SKE = skewness(datinroi);
    xlsout{i+1,1} = files(i).name(1:end-4);
    xlsout{i+1,2} = voxnum;
    xlsout{i+1,3} = vols;
    xlsout{i+1,4} = minv;
    xlsout{i+1,5} = maxv;
    xlsout{i+1,6} = meanv;
    xlsout{i+1,7} = stdv;
    xlsout{i+1,8} = medianv;
    xlsout{i+1,9} = Q1v;
    xlsout{i+1,10} = Q3v;
    xlsout{i+1,11} = KUR;
    xlsout{i+1,12} = SKE;
%     xlsout(i+1,2) = voxnum;
%     xlsout(i+1,3) = vols;
%     xlsout(i+1,4) = minv;
%     xlsout(i+1,5) = maxv;
%     xlsout(i+1,6) = meanv;
%     xlsout(i+1,7) = stdv;
%     xlsout(i+1,8) = medianv;
%     xlsout(i+1,9) = Q1v;
%     xlsout(i+1,10) = Q3v;
%     xlsout(i+1,11) = KUR;
%     xlsout(i+1,12) = SKE;
    save([OutHistdir,filesep,files(i).name(1:end-4),'.mat'],'datinroi');
end
save([OutHistdir,filesep,'HistgramInfor.mat'],'xlsout');
try
    xlswrite([OutHistdir,filesep,'HistgramInfor.xls'],xlsout);
catch
    
end
   