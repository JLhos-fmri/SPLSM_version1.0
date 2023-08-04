function ClincFreq_D2V(INDIR,OutD2Vdir,voxelsize,Ventricledir,Ventriclematdir)
[Vven,Dven] = Dynamic_read_dir_NIFTI(Ventricledir);
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;
load(Ventriclematdir);
Dvens = reshape(Dven,Vven.dim);
for i = 1:size(ROItab,1)
    indexs = find(Dvens==ROItab{i,1});
    [ix iy iz] = ind2sub(Vven.dim,indexs);
    ind{i,1} = [ix iy iz];
    clear indexs ix iy iz;
end
% inum = 1;
for i = 1:size(ROItab,1)
    xlsinfo{1,i+1} = ROItab{i,2};
end
xlsinfoA = xlsinfo;
xlsinfomm = xlsinfo;
xlsinfommA = xlsinfo;
for isub = 1:length(namelist)
    Du = reshape(full(D(:,isub)),V(1).dim);
    Coninf = bwconncomp(Du);
    nums = Coninf.NumObjects;
    if nums>0
        for i = 1:nums
            Dattemp = zeros(size(Du));
            Dattemp(Coninf.PixelIdxList{i}) = 1;
            stats = regionprops3(Dattemp);
            Cent(i,:) = stats.Centroid;
            %             Vol(i,:) = stats.Volume;
        end
        CentOut = [round(Cent(:,2)),round(Cent(:,1)),round(Cent(:,3))];
        for i = 1:nums
            for j = 1:size(ROItab,1)
                inds = ind{j,1};
                dist(i,j) = min(sqrt((CentOut(i,1)-inds(:,1)).^2+(CentOut(i,2)-inds(:,2)).^2+(CentOut(i,3)-inds(:,3)).^2));
            end
        end
        for i = 1:nums
            [ixD iyD izD] = ind2sub(Vven.dim,Coninf.PixelIdxList{i});
            for j = 1:size(ROItab,1)
                parfor k = 1:length(ixD)
                    inds = ind{j,1};
                    distV(k,j) = min(sqrt((ixD(k)-inds(:,1)).^2+(iyD(k)-inds(:,2)).^2+(izD(k)-inds(:,3)).^2));
                end
            end
            distVox(i,:) = min(distV,[],1);
            clear distV
        end
        
        prename = namelist{isub};
        save([OutD2Vdir,filesep,prename(1:end-4),'.mat'],'dist','distVox','voxelsize')
        xlsinfo{isub+1,1} = prename;
        xlsinfoA{isub+1,1} = prename;
        for i = 1:nums
            for j = 1:size(ROItab,1)
                xlsinfo{isub+1,(i-1)*size(ROItab,1)+1+j} = dist(i,j);
                xlsinfoA{isub+1,(i-1)*size(ROItab,1)+1+j} = distVox(i,j);
            end
        end
        
        xlsinfomm{isub+1,1} = prename;
        xlsinfommA{isub+1,1} = prename;
        for i = 1:nums
            for j = 1:size(ROItab,1)
                xlsinfomm{isub+1,(i-1)*size(ROItab,1)+1+j} = dist(i,j)*voxelsize;
                xlsinfommA{isub+1,(i-1)*size(ROItab,1)+1+j} = distVox(i,j)*voxelsize;
            end
        end
        clear dist distVox
    else
        prename = namelist{isub};
        xlsinfo{isub+1,1} = prename;
        xlsinfomm{isub+1,1} = prename;
        xlsinfoA{isub+1,1} = prename;
        xlsinfommA{isub+1,1} = prename;
    end
end
save([OutD2Vdir,filesep,'Dist2Ven.mat'],'xlsinfo','xlsinfomm','xlsinfoA','xlsinfommA')
try
    xlswrite([OutD2Vdir,filesep,'Dist2Ven_voxel.xls'],xlsinfo)
    xlswrite([OutD2Vdir,filesep,'Dist2Ven_mm.xls'],xlsinfomm)
    xlswrite([OutD2Vdir,filesep,'Dist2VenAll_voxel.xls'],xlsinfoA)
    xlswrite([OutD2Vdir,filesep,'Dist2VenAll_mm.xls'],xlsinfommA)
catch
    
end
end
