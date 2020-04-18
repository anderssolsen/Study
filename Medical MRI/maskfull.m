function mask = maskfull(vol,threshold)


mask = ones(size(vol));

%remove from 160 up

for slice = 1:240

data = squeeze(vol(slice,:,:));
for i=1:240
    
    j = 1;
    while j
        if j>160
            j = [];
        elseif data(i,j)<threshold
            mask(slice,i,j)=0;
            j = j+1;
        else, mask(slice,i,j:j+12)=0;
            j = [];
%             mask(i,j-3:j-1)=1;
        end
    end
    
    k = 1;
    while k
        if k>160
            k = [];
        elseif data(i,161-k)<threshold
            mask(slice,i,161-k)=0;
            k = k+1;
        else, mask(slice,i,161-k-13:161-k)=0;
            k = [];
        end
    end
    
    
end
end
mask=mask(1:240,1:240,1:160);
mask(160:end,:,:) = 0;
mask = logical(mask);
end