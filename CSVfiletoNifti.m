function CSVfiletoNifti(dirtempschi001,Vatl,Datl)
Datl(isnan(Datl)) = 0;
indfordatout = unique(Datl);
indfordatout(indfordatout==0) = [];
csvfiles = dir([dirtempschi001,filesep,'*.csv']);
for i = 1:length(csvfiles)
    outname = [dirtempschi001,filesep,csvfiles(i).name(1:end-4),'.nii'];
    P05 = csvread([dirtempschi001,filesep,csvfiles(i).name]);
    DatOut = zeros(Vatl(1).dim);
    for iout = 1:length(indfordatout)
        DatOut(Datl==indfordatout(iout)) = P05(iout);
    end
    DynamicBC_write_NIFTI(DatOut,Vatl(1),outname)
end
end