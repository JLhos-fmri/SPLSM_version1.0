function CFA_MakeRes_Clustercorrected(C05,C05P,C05N,PmapDir,THRESHOLD,VQ,Voxdir2out05)

[VP DP] = Dynamic_read_dir_NIFTI(PmapDir);
[pth nam ext] = fileparts(PmapDir);

THR1P = find(DP>1-THRESHOLD);
THR1N = find(DP<THRESHOLD&DP>0);
THR1 = zeros(size(DP));
THR1_P = zeros(size(DP));
THR1_N = zeros(size(DP));
THR1(THR1P) = DP(THR1P);
THR1(THR1N) = DP(THR1N)-1;
THR1_P(THR1P) = DP(THR1P);
THR1_N(THR1N) = DP(THR1N)-1;

THR1_Reshape = reshape(THR1,VQ.dim);
THR1_P_Reshape = reshape(THR1_P,VQ.dim);
THR1_N_Reshape = reshape(THR1_N,VQ.dim);

Coninf = bwconncomp(logical(THR1_Reshape));
nums = Coninf.NumObjects;
P05num = 0;
P01num = 0;
P005num = 0;
P001num = 0;
if nums>0
    for i = 1:nums
        clusnum = length(Coninf.PixelIdxList{i});
        if clusnum>=C05.P0001
            P001num = P001num+1;
            THR1_C001{P001num} = Coninf.PixelIdxList{i};
            P005num = P005num+1;
            THR1_C005{P005num} = Coninf.PixelIdxList{i};
            P01num = P01num+1;
            THR1_C01{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05.P0005
            P005num = P005num+1;
            THR1_C005{P005num} = Coninf.PixelIdxList{i};
            P01num = P01num+1;
            THR1_C01{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05.P001
            P01num = P01num+1;
            THR1_C01{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05.P005
            P05num = P05num+1;
            THR1_C05{P05num} = Coninf.PixelIdxList{i};
        end
    end
end
if P001num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P001num)
        OutMap(THR1_C001{i}) = THR1_Reshape(THR1_C001{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Both_C0001.nii'])
end
if P005num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P005num)
        OutMap(THR1_C005{i}) = THR1_Reshape(THR1_C005{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Both_C0005.nii'])    
end
if P01num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P01num)
        OutMap(THR1_C01{i}) = THR1_Reshape(THR1_C01{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Both_C001.nii'])    
end
if P05num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P05num)
        OutMap(THR1_C05{i}) = THR1_Reshape(THR1_C05{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Both_C005.nii'])    
end

%% pos

Coninf = bwconncomp(logical(THR1_P_Reshape));
nums = Coninf.NumObjects;
P05num = 0;
P01num = 0;
P005num = 0;
P001num = 0;
if nums>0
    for i = 1:nums
        clusnum = length(Coninf.PixelIdxList{i});
        if clusnum>=C05P.P0001
            P001num = P001num+1;
            THR1_C001P{P001num} = Coninf.PixelIdxList{i};
            P005num = P005num+1;
            THR1_C005P{P005num} = Coninf.PixelIdxList{i};
            P01num = P01num+1;
            THR1_C01P{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05P{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05P.P0005
            P005num = P005num+1;
            THR1_C005P{P005num} = Coninf.PixelIdxList{i};
            P01num = P01num+1;
            THR1_C01P{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05P{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05P.P001
            P01num = P01num+1;
            THR1_C01P{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05P{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05P.P005
            P05num = P05num+1;
            THR1_C05P{P05num} = Coninf.PixelIdxList{i};
        end
    end
end
if P001num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P001num)
        OutMap(THR1_C001P{i}) = THR1_Reshape(THR1_C001P{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Pos_C0001.nii'])
end
if P005num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P005num)
        OutMap(THR1_C005P{i}) = THR1_Reshape(THR1_C005P{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Pos_C0005.nii'])    
end
if P01num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P01num)
        OutMap(THR1_C01P{i}) = THR1_Reshape(THR1_C01P{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Pos_C001.nii'])    
end
if P05num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P05num)
        OutMap(THR1_C05P{i}) = THR1_Reshape(THR1_C05P{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Pos_C005.nii'])    
end

%% Neg

Coninf = bwconncomp(logical(THR1_N_Reshape));
nums = Coninf.NumObjects;
P05num = 0;
P01num = 0;
P005num = 0;
P001num = 0;
if nums>0
    for i = 1:nums
        clusnum = length(Coninf.PixelIdxList{i});
        if clusnum>=C05N.P0001
            P001num = P001num+1;
            THR1_C001N{P001num} = Coninf.PixelIdxList{i};
            P005num = P005num+1;
            THR1_C005N{P005num} = Coninf.PixelIdxList{i};
            P01num = P01num+1;
            THR1_C01N{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05N{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05N.P0005
            P005num = P005num+1;
            THR1_C005N{P005num} = Coninf.PixelIdxList{i};
            P01num = P01num+1;
            THR1_C01N{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05N{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05N.P001
            P01num = P01num+1;
            THR1_C01N{P01num} = Coninf.PixelIdxList{i};
            P05num = P05num+1;
            THR1_C05N{P05num} = Coninf.PixelIdxList{i};
        elseif clusnum>=C05N.P005
            P05num = P05num+1;
            THR1_C05N{P05num} = Coninf.PixelIdxList{i};
        end
    end
end
if P001num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P001num)
        OutMap(THR1_C001N{i}) = THR1_Reshape(THR1_C001N{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Neg_C0001.nii'])
end
if P005num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P005num)
        OutMap(THR1_C005N{i}) = THR1_Reshape(THR1_C005N{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Neg_C0005.nii'])    
end
if P01num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P01num)
        OutMap(THR1_C01N{i}) = THR1_Reshape(THR1_C01N{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Neg_C001.nii'])    
end
if P05num>0
    OutMap = zeros(VQ.dim);
    for i = 1:length(P05num)
        OutMap(THR1_C05N{i}) = THR1_Reshape(THR1_C05N{i});        
    end
    DynamicBC_write_NIFTI(OutMap,VQ,[Voxdir2out05,nam,'marked_Neg_C005.nii'])    
end

% DynamicBC_write_NIFTI(reshape(THR1,VQ.dim),VQ,[Voxdir2out05,nam,'marked_P.nii']);

end