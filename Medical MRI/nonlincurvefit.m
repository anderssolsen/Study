function [PD,T1] = nonlincurvefit(slices,mask,alpha,TR)
% This function finds PD and T1 values voxel-wise on an array of minimum 3
% slices, either coronal, saggital or transversal

% preprocessing
slices  = squeeze(slices);
szorig  = size(slices);
Nvoxels = szorig(1)*szorig(2);
slices  = reshape(slices,Nvoxels,szorig(3));

mask    = reshape(mask,Nvoxels,1);
Nmask   = sum(mask);

slices  = double(slices(mask,:)); %%% Use only voxels within mask
slices  = slices./max(slices(:)); %normalize to 0 and 1


%%% Curvefit
x = zeros(Nmask,2); %preallocation

fun = @(x,alpha)sind(alpha).*x(1).*(1-exp(-TR./x(2)))./(1-cosd(alpha).*exp(-TR./x(2)));
x0 = [0.1,0.1];
opts = optimset('Display','off');
lowbound = [0,0];

tic
for a = 1:Nmask
    x(a,:) = lsqcurvefit(fun,x0,alpha,slices(a,:),lowbound,[],opts);

    if ismember(a,[10000:10000:Nmask])
        disp(['Done with nr ',num2str(a),' out of ',num2str(Nmask)])
        toc
    end
end

%%% Output
PD = zeros(Nvoxels,1);
PD(mask) = x(:,1);
T1 = zeros(Nvoxels,1);
T1(mask) = x(:,2);

PD = reshape(PD,szorig(1),szorig(2));
T1 = reshape(T1,szorig(1),szorig(2));

end
