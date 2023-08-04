function [PvalpermV005C] = getclusterinfo2(PvalpermV005)


Coninf = bwconncomp(PvalpermV005);
nums = Coninf.NumObjects;
if nums>0
    for i = 1:nums
        PvalpermV005C(i,1) = length(Coninf.PixelIdxList{i});
    end
else
    PvalpermV005C = 0;
end


end