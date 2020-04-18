%%% Initialize
clear
directory = 'C:\Users\ander\Documents\10. Semester\Medical MRI\ExerciseData\ass1';
datadir = strcat(directory,'\Data\');
addpath(genpath(directory))
cd(directory)


%% load data and extract information
i1 = get_mri_data(datadir,{'X98279_WIP_sT1W_3D_FFE_31547_4_1'},1);
i2 = get_mri_data(datadir,{'X98279_WIP_sT1W_3D_FFE_31547_5_1'},1);
i3 = get_mri_data(datadir,{'X98279_WIP_sT1W_3D_FFE_31547_6_1'},1);
i4 = get_mri_data(datadir,{'X98279_WIP_sT1W_3D_FFE_31547_7_1'},1);
alpha = [i1.info.TipAngle,i2.info.TipAngle,i3.info.TipAngle,i4.info.TipAngle];
nVoxels = size(i1.data,1)*size(i1.data,2)*size(i1.data,3);
TR = i1.info.RepetitionTime;
numims = 4;

%% Decide which slices to use
slice = 90;
imagesc([i1.data(:,:,slice), i2.data(:,:,slice), ...
    i3.data(:,:,slice), i4.data(:,:,slice)]);colormap gray;axis image

data = cat(4,i1.data,i2.data,i3.data,i4.data);
dataM = fliplr(data(:,:,80,:));
montage ( dataM ,'Size',[1,4],'DisplayRange',[0 1600])

%% perform (simple) masking on selected slices
sag80 = i2.data(:,:,80);
cor80 = i2.data(80,:,:);
trans140 = i2.data(:,140,:);

threshold = 150;
masksag80 = masksimple(sag80,threshold);
maskcor80 = masksimple(cor80,threshold);
masktrans140 = masksimple(trans140,threshold);

figure;
subplot(1,2,1)
imagesc(sag80);colormap gray
subplot(1,2,2)
imagesc(masksag80);colormap gray

%% perform PD and T1 eval on the saggital slice

sag80 = data(:,:,90,:);

[PDsag,T1sag] = nonlincurvefit(sag80,masksag80,alpha,TR);

%%% Remove values above a threshold
PDsag2 = PDsag;
PDsag2(PDsag2>20)=0;
T1sag2 = T1sag;
T1sag2(T1sag2>6)=0;
figure;
subplot(1,2,1);imagesc(PDsag2);colormap gray; axis image
subplot(1,2,2);imagesc(T1sag2);colormap gray; axis image

%%%%%%% Normalizing PD using polygon mask on ventricle (saved)
load polycsf.mat %polygon for ventricle
ROIcsf = createMask(polycsf);
PDcsf = mean(PDsag(ROIcsf));
T1csf = mean(T1sag(ROIcsf));
PDsag3 = PDsag2./PDcsf;

%%%% find corpus callosum mask values
load polyWM.mat
ROIWM = createMask(polyWM);
PDWM = mean(PDsag3(ROIWM));
T1WM = mean(T1sag2(ROIWM));

%% Do the same on the other sides

cor80 = data(90,:,:,:);
trans140 = data(:,140,:,:);
[PDcor,T1cor] = nonlincurvefit(cor80,maskcor80,alpha,TR);
[PDtrans,T1trans] = nonlincurvefit(trans140,masktrans140,alpha,TR);

PDcor2 = PDcor./PDcsf;
PDcor2(PDcor2>20)=0;
T1cor2 = T1cor;
T1cor2(T1cor2>6)=0;

PDtrans2 = PDtrans./PDcsf;
PDtrans2(PDtrans2>20)=0;
T1trans2 = T1trans;
T1trans2(T1trans2>6)=0;


PDdata = [fliplr(PDsag3),flipud(PDcor2),PDtrans2];
figure,imagesc(PDdata,[0,3]);colormap gray;axis image;axis off
T1data = [fliplr(T1sag2),flipud(T1cor2),T1trans2];
figure,imagesc(T1data,[0,6]);colormap gray;axis image;axis off

figure,imagesc(flipud(T1cor2));colormap gray;axis image

%%% Find values within putamen
load polyGM.mat
ROIGM = createMask(polyGM);
PDGM = mean(PDcor2(ROIGM))
T1GM = mean(T1cor2(ROIGM))
%% Do on full volume

mask = maskfull(i2.data,150); %approved
PDfull = zeros(size(i2.data));
T1full = zeros(size(i2.data));

for slice = 1:160
    dat = squeeze(data(slice,:,:,:));
    mas = squeeze(mask(slice,:,:));
    [PDfull(slice,:,:),T1full(slice,:,:)] = nonlincurvefit(dat,mas,alpha,TR);
%changed disp in fun to 10000
disp(['Done with slice ',num2str(slice),' out of 160'])
end


PDfull2 = PDfull;
PDfull2(PDfull2>10)=0;
T1full2 = T1full;
T1full2(T1full2>6)=0;

slicei(PDfull2);colormap gray;
slicei(T1full2);colormap gray;

figure;imagesc(squeeze(PDfull2(:,:,80)));axis image;colormap gray

slice90 = flipud(squeeze(T1full2(90,:,:)));
ROIGM = createMask(polyGM);
T1GM = mean(slice90(ROIGM));




