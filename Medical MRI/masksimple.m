function mask = masksimple(slice,threshold)
% Construct a simple mask by going from left to right and right to left in
% an image and output 0 if the voxel value is less than threshold,
% otherwise zero. For every line, the algorithm stops upon first encounter
% of a voxel value greater than threshold, to avoid masking out the inside
% of the head.
% To achieve good masking, select an image (i3) with high voxel values at
% the skin, and a high threshold to mask out everything outside the head.
slice = squeeze(slice);
topdown = size(slice,1);
leftright = size(slice,2);

mask = ones(size(slice));

for i=1:topdown
    
    j = 1;
    while j
        if j>leftright
            j = [];
        elseif slice(i,j)<threshold
            mask(i,j)=0;
            j = j+1;
        else
            j = [];
        end
    end
    
    k = 1;
    while k
        if k>leftright
            k = [];
        elseif slice(i,leftright+1-k)<threshold
            mask(i,leftright+1-k)=0;
            k = k+1;
        else
            k = [];
        end
    end
    
    
end
mask = logical(mask);
end