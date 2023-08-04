function CFA_CheckData
% uiwait(msgbox({'Currently, it only used for two conditions:',...
%     '1. Data after regrouped but before normalization.',...
%     '2. Data after normalization.'}))
PG = uigetdir(pwd,'Check Data');
% PGout = uigetdir(pwd,'Result Picture directory');
[pth nam ext] = fileparts(PG);
PGout = [pth,filesep,'ROICheckBeforeNormalization'];
mkdir(PGout);
ansct = questdlg('CT mode','CT mode','yes','no','no');
sfold = dir(PG);
cols = [gray;[ones(64,1),zeros(64,1),zeros(64,1)]];
for i = 1:length(sfold)-2
    dirtemps = [PG,filesep,sfold(i+2).name,filesep];
    fileimage = [dirtemps,'image.nii'];
    fileroi = [dirtemps,'ROI.nii'];
    fileT1 = [dirtemps,'T1.nii'];
    fileother = [dirtemps,'OtherMode.nii'];
    try
        if ~isempty(dir(fileimage)) %
            [VI DI] = Dynamic_read_dir_NIFTI(fileimage);
            [VR DR] = Dynamic_read_dir_NIFTI(fileroi);
            if any(sum(VI.mat-VR.mat))
                disp([sfold(i+2),': ROI not match image']);
                continue;
            end
            DI(isnan(DI)) = 0;
            DI(isinf(DI)) = 0;
            DR(isnan(DR)) = 0;
            DR(isinf(DR)) = 0;
            
            Dimage = reshape(DI,VI.dim);
            Droi = reshape(DR,VR.dim);
            Dimgrang = unique(DI(DI>0));
            %             Dimgrang = DI(DI>0);
            if strcmp(ansct,'no')
                [ix iy] = sort(Dimgrang);
                Drangimage(1) = ix(round(length(ix)*0.05));
                Drangimage(2) = ix(round(length(ix)*0.95));
            else
                Drangimage(1) = 10;
                Drangimage(2) = 100;
            end
            Dshow = Dimage;
            Dshow(Dimage<Drangimage(1)) = Drangimage(1);
            Dshow(Dimage>Drangimage(2)) = Drangimage(2);
            Dshow(Dshow>=Drangimage(1)) = (Dshow(Dshow>=Drangimage(1))-Drangimage(1))/(Drangimage(2)-Drangimage(1))*0.95+0.025;
            Dshow2 = reshape(Dshow,VI.dim);
            
            Dshow2(Droi>0) = 1.5;
            TotalSlice = VR.dim(3);
            if TotalSlice>72
                Coninf = bwconncomp(Droi);
                for iroi = 1:Coninf.NumObjects
                    pixlist = Coninf.PixelIdxList{iroi};
                    dattemp = zeros(VR.dim);
                    dattemp(pixlist) = 1;
                    stats = regionprops3(dattemp);
                    Cent = stats.Centroid;
                    Pos(iroi,1:3) = round(Cent([2,1,3]));
                    if round(Cent(3))<1
                        importslice(iroi,1) = 1;
                    elseif round(Cent(3))>=VR.dim(3)
                        importslice(iroi,1) = VR.dim(3);
                    else
                        importslice(iroi,1) = round(Cent(3));
                    end
                end
                UsedSlice = unique(importslice);
                Sliceshow = floor(TotalSlice/72);
                Sliceadd = floor((TotalSlice-72*Sliceshow)/2);
                SliceShowind = Sliceadd+1:Sliceshow:(TotalSlice-Sliceadd);
                
                [datint,indse,indse2] = intersect(UsedSlice,SliceShowind);
                indc = UsedSlice;
                indc(indse) = [];
                indraw = SliceShowind;
                %             indraw = [];
                %             SliceShowfinal = SliceShowind(indse2);
                if ~isempty(indc)
                    for inum = 1:length(indc)
                        if length(indraw)>0
                            [minv mindis] = min(abs(indc(inum)-indraw));
                            indraw(mindis) = [];
                            indraw(end+1) = indc(inum);
                            %                         SliceShowfinal(end+1) = indc(inum);
                        end
                    end
                end
                Sliceshowix = sort(indraw);
                %             Sliceshowix = sort(SliceShowfinal);
                hfig = figure('unit','norm',...
                    'pos',[0.1,0.1,0.8,0.8],...
                    'name',sfold(i+2).name);
                for inum = 1:min(length(Sliceshowix),72)
                    subplot(8,9,inum);imagesc(rot90(squeeze(Dshow2(:,:,Sliceshowix(inum))),3),[0,2]);
                    colormap(cols);
                    axis off;
                end
            else
                hfig = figure('unit','norm',...
                    'pos',[0.1,0.1,0.8,0.8],...
                    'name',sfold(i+2).name);
                if TotalSlice<=20
                    for islice = 1:TotalSlice
                        subplot(4,5,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        colormap(cols);
                        axis off
                    end
                elseif TotalSlice<=30
                    for islice = 1:TotalSlice
                        subplot(5,6,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        colormap(cols);
                        axis off
                    end
                elseif TotalSlice<=42
                    for islice = 1:TotalSlice
                        subplot(6,7,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        colormap(cols);
                        axis off
                    end
                elseif TotalSlice<=56
                    for islice = 1:TotalSlice
                        subplot(7,8,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        colormap(cols);
                        axis off
                    end
                elseif TotalSlice<=72
                    for islice = 1:TotalSlice
                        subplot(8,9,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        colormap(cols);
                        axis off
                    end
                end
            end
        else % with T1 mode
            [vT1, dT1] = Dynamic_read_dir_NIFTI(fileT1);
            [VI DI] = Dynamic_read_dir_NIFTI(fileother);
            [VR DR] = Dynamic_read_dir_NIFTI(fileroi);
            
            if any(sum(VI.mat-VR.mat))
                disp([sfold(i+2),': ROI not match image']);
                continue;
            end
            DI(isnan(DI)) = 0;
            DI(isinf(DI)) = 0;
            DR(isnan(DR)) = 0;
            DR(isinf(DR)) = 0;
            
            Dimage = reshape(DI,VI.dim);
            Droi = reshape(DR,VR.dim);
            Dimgrang = unique(DI(DI>0));
            %             Dimgrang = DI(DI>0);
            [ix iy] = sort(Dimgrang);
            %         Drangimage(1) = ix(round(length(ix)*0.05));
            %         Drangimage(2) = ix(round(length(ix)*0.95));
            %         Dshow = Dimage;
            %         Dshow(Dimage<Drangimage(1)) = Drangimage(1);
            %         Dshow(Dimage>Drangimage(2)) = Drangimage(2);
            %         Dshow2 = reshape(Dshow,VI.dim);
            
            
            Drangimage(1) = ix(round(length(ix)*0.05));
            Drangimage(2) = ix(round(length(ix)*0.95));
            Dshow = Dimage;
            Dshow(Dimage<Drangimage(1)) = Drangimage(1);
            Dshow(Dimage>Drangimage(2)) = Drangimage(2);
            Dshow(Dshow>=Drangimage(1)) = (Dshow(Dshow>=Drangimage(1))-Drangimage(1))/(Drangimage(2)-Drangimage(1))*0.95+0.025;
            Dshow2 = reshape(Dshow,VI.dim);
            
            Dshow2(Droi>0) = 1.5;
            
            TotalSlice = VR.dim(3);
            if TotalSlice>72
                Coninf = bwconncomp(Droi);
                for iroi = 1:Coninf.NumObjects
                    pixlist = Coninf.PixelIdxList{iroi};
                    dattemp = zeros(VR.dim);
                    dattemp(pixlist) = 1;
                    stats = regionprops3(dattemp);
                    Cent = stats.Centroid;
                    Pos(iroi,1:3) = round(Cent([2,1,3]));
                    if round(Cent(3))<1
                        importslice(iroi,1) = 1;
                    elseif round(Cent(3))>=VR.dim(3)
                        importslice(iroi,1) = VR.dim(3);
                    else
                        importslice(iroi,1) = round(Cent(3));
                    end
                end
                UsedSlice = unique(importslice);
                Sliceshow = floor(TotalSlice/72);
                Sliceadd = floor((TotalSlice-72*Sliceshow)/2);
                SliceShowind = Sliceadd+1:Sliceshow:(TotalSlice-Sliceadd);
                
                [datint,indse,indse2] = intersect(UsedSlice,SliceShowind);
                indc = UsedSlice;
                indc(indse) = [];
                indraw = SliceShowind;
                indraw(indse2) = [];
                SliceShowfinal = SliceShowind(indse2);
                
                
                %             SliceShowfinal = SliceShowind(indse2);
                if ~isempty(indc)
                    for inum = 1:length(indc)
                        if length(indraw)>0
                            [minv mindis] = min(abs(indc(inum)-indraw));
                            indraw(mindis) = [];
                            indraw(end+1) = indc(inum);
                            %                         SliceShowfinal(end+1) = indc(inum);
                        end
                    end
                end
                Sliceshowix = sort(indraw);
                %             Sliceshowix = sort(SliceShowfinal);
                
                hfig = figure('unit','norm',...
                    'pos',[0.1,0.1,0.8,0.8],...
                    'name',sfold(i+2).name);
                for inum = 1:length(Sliceshowix)
                    subplot(8,9,inum);imagesc(rot90(squeeze(Dshow2(:,:,Sliceshowix(inum))),3),[0,2]);
                    axis off
                end
            else
                hfig = figure('unit','norm',...
                    'pos',[0.1,0.1,0.8,0.8],...
                    'name',sfold(i+2).name);
                if TotalSlice<=20
                    for islice = 1:TotalSlice
                        subplot(4,5,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        axis off
                    end
                elseif TotalSlice<=30
                    for islice = 1:TotalSlice
                        subplot(5,6,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        axis off
                    end
                elseif TotalSlice<=42
                    for islice = 1:TotalSlice
                        subplot(6,7,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        axis off
                    end
                elseif TotalSlice<=56
                    for islice = 1:TotalSlice
                        subplot(7,8,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        axis off
                    end
                elseif TotalSlice<=72
                    for islice = 1:TotalSlice
                        subplot(8,9,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                        axis off
                    end
                end
            end
        end
        colormaps = [gray;ones(64,1),zeros(64,1),zeros(64,1)];
        colormap(colormaps)
        print(hfig,'-djpeg','-r300',[PGout,filesep,sfold(i+2).name,'.jpg']);
        close(hfig)
    catch
        disp(['Wrong ROI in ',sfold(i+2).name])
    end
end
%%
end