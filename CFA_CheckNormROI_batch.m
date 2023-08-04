function CFA_CheckNormROI_batch(PG,CTmode)
% PG = uigetdir(pwd,'Check Data(after normalize)');
dir1 = [PG,filesep,'NormBG',filesep];
dir2 = [PG,filesep,'NormROI',filesep];
dir3 = [PG,filesep,'NormBGT1',filesep];
file1 = dir([dir1,'*.nii']);
file2 = dir([dir2,'*.nii']);
file3 = dir([dir3,'*.nii']);

% ansct = questdlg('CT mode','CT mode','yes','no','no');

PGout = [PG,filesep,'CheckNormPicture',filesep];
mkdir(PGout)
colos = [gray(64);parula(64);[ones(64,1),zeros(64,1),zeros(64,1)]];
for i = 1:length(file1)
    [v1 d1] = Dynamic_read_dir_NIFTI([dir1,file1(i).name]);
    d1(isnan(d1)) = 0;
    d1(isinf(d1)) = 0;
    %%
    %     Dimgrang = unique(DI(DI>0));
    %     [ix iy] = sort(Dimgrang);
    %     Drangimage(1) = ix(round(length(ix)*0.05));
    %     Drangimage(2) = ix(round(length(ix)*0.95));
    %%
    Drang = unique(d1(d1>0));
    [ix iy] = sort(Drang);
    if ~CTmode
        dx1 = ix(round(length(ix)*0.05));
        dx2 = ix(round(length(ix)*0.95));
    else
        dx1 = 1;
        dx2 = 100;
    end
    d1(d1<dx1&d1>0) = dx1;
    d1(d1>dx2) = dx2;
    d1(d1>=dx1) = (d1(d1>=dx1)-dx1)/(dx2-dx1)*0.95+0.025;
    DAT1 = reshape(d1,v1.dim);
    TotalSlice = size(DAT1,3);
    
    [v2 d2] = Dynamic_read_dir_NIFTI([dir2,file2(i).name]);
    d2(isnan(d2)) = 0;
    d2(isinf(d2)) = 0;
%     d2(d2>=1) = 2.5;
    DAT2 = reshape(d2,v2.dim);
    DAT1(DAT2>=0.5) = 2.5;
    if ~isempty(file3)
        [v3 d3] = Dynamic_read_dir_NIFTI([dir3,file3(i).name]);
        d3(isnan(d3)) = 0;
        d3(isinf(d3)) = 0;
        Drang = unique(d3(d3>0));
        [ix iy] = sort(Drang);
        dx1 = ix(round(length(ix)*0.05));
        dx2 = ix(round(length(ix)*0.95));
        d3(d3<dx1&d3>0) = dx1;
        d3(d3>dx2) = dx2;
        d3(d3>=dx1) = (d3(d3>=dx1)-dx1)/(dx2-dx1)*0.95+1.025;
        DAT3 = reshape(d3,v3.dim);
        %%
        %         Dimgrang = unique(DI(DI>0));
        %         [ix iy] = sort(Dimgrang);
        %         Drangimage(1) = ix(round(length(ix)*0.05));
        %         Drangimage(2) = ix(round(length(ix)*0.95));
        %%
    end
    Dshow2 = DAT1;
    if TotalSlice>72
        Coninf = bwconncomp(DAT2);
        for iroi = 1:Coninf.NumObjects
            pixlist = Coninf.PixelIdxList{iroi};
            dattemp = zeros(v1.dim);
            dattemp(pixlist) = 1;
            stats = regionprops3(dattemp);
            Cent = stats.Centroid;
            Pos(iroi,1:3) = round(Cent([2,1,3]));
            if round(Cent(3))<1
                importslice(iroi,1) = 1;
            elseif round(Cent(3))>=v1.dim(3)
                importslice(iroi,1) = v1.dim(3);
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
        %         SliceShowfinal = indraw;
        %         SliceShowfinal = SliceShowind(indse2);
        if ~isempty(indc)
            for inum = 1:length(indc)
                if length(indraw)>0
                    [minv mindis] = min(abs(indc(inum)-indraw));
                    indraw(mindis) = [];
                    indraw(end+1) = indc(inum);
                    %                     SliceShowfinal(end+1) = indc(inum);
                end
            end
        end
        Sliceshowix = sort(indraw);
        %         Sliceshowix = sort(SliceShowfinal);
        hfig = figure('unit','norm',...
            'pos',[0.1,0.1,0.8,0.8],...
            'name',file1(i).name);
        for inum = 1:min(length(Sliceshowix),72)
            subplot(8,9,inum);imagesc(rot90(squeeze(Dshow2(:,:,Sliceshowix(inum))),1),[0,3]);
            %             subplot(8,9,inum);imagesc(rot90(squeeze(Dshow2(:,:,Sliceshowix(inum))),3),[0,2]);
            hold on;
            if ~isempty(file3)
                h1 = imagesc(rot90(DAT3(:,:,Sliceshowix(inum)),1),[0,3]);
                %             h1 = imagesc(rot90(DAT3(:,:,Sliceshowix(inum)),3),[0,3]);
                set(h1,'AlphaData',0.5);
            end
            colormap(colos);
            axis off;
        end
    else
        hfig = figure('unit','norm',...
            'pos',[0.1,0.1,0.8,0.8],...
            'name',file1(i).name);
        if TotalSlice<=20
            for islice = 1:TotalSlice
                subplot(4,5,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),1),[0,3]);
                %                 subplot(4,5,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                hold on;
                if ~isempty(file3)
                    h1 = imagesc(rot90(DAT3(:,:,islice),1),[0,3]);
                    %                 h1 = imagesc(rot90(DAT3(:,:,islice),3),[0,3]);
                    set(h1,'AlphaData',0.5);
                end
                colormap(colos);
                axis off;
            end
        elseif TotalSlice<=30
            for islice = 1:TotalSlice
                subplot(5,6,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),1),[0,3]);
                %                 subplot(5,6,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                hold on;
                if ~isempty(file3)
                    h1 = imagesc(rot90(DAT3(:,:,islice),1),[0,3]);
                    %                 h1 = imagesc(rot90(DAT3(:,:,islice),3),[0,3]);
                    set(h1,'AlphaData',0.5);
                end
                colormap(colos);
                axis off;
            end
        elseif TotalSlice<=42
            for islice = 1:TotalSlice
                subplot(6,7,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),1),[0,3]);
                %                 subplot(6,7,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                hold on;
                if ~isempty(file3)
                    h1 = imagesc(rot90(DAT3(:,:,islice),1),[0,3]);
                    %                 h1 = imagesc(rot90(DAT3(:,:,islice),3),[0,3]);
                    set(h1,'AlphaData',0.5);
                end
                colormap(colos);
                axis off;
            end
        elseif TotalSlice<=56
            for islice = 1:TotalSlice
                subplot(7,8,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),1),[0,3]);
                %                 subplot(7,8,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                hold on;
                if ~isempty(file3)
                    h1 = imagesc(rot90(DAT3(:,:,islice),1),[0,3]);
                    %                 h1 = imagesc(rot90(DAT3(:,:,islice),3),[0,3]);
                    set(h1,'AlphaData',0.5);
                end
                colormap(colos);
                axis off;
            end
        elseif TotalSlice<=72
            for islice = 1:TotalSlice
                subplot(8,9,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),1),[0,3]);
                %                 subplot(8,9,islice);imagesc(rot90(squeeze(Dshow2(:,:,islice)),3),[0,2]);
                hold on;
                if ~isempty(file3)
                    h1 = imagesc(rot90(DAT3(:,:,islice),1),[0,3]);
                    %                 h1 = imagesc(rot90(DAT3(:,:,islice),3),[0,3]);
                    set(h1,'AlphaData',0.5);
                end
                colormap(colos);
                axis off;
            end
        end
    end
    print(hfig,'-djpeg','-r300',[PGout,file1(i).name(1:end-4),'.jpg']);
    close(hfig);
    
end
disp('check Normalized ROIs finished!')
end