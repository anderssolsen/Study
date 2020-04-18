function hf=slicei(im,xyzt)
% Usage: hf=slicei(im,[xyzt]);
% Given "im", a three-dimensional field that may also vary over time,
% slicei visualizes the data via an interactive figure with orthographic
% projections and a three-dimensional slice plot. The user can move around
% in all three directions by dragging. If the field varies over time, the
% user can also change the time, play a movie of the data, and save a movie
% of the data to disk. If "im" is a scalar, provide it as a three-
% dimensional rectangular array with data arranged for a grid produced by
% meshgrid.m (not ndgrid.m). If "im" is a vector, provide it as a cell
% array of three three-dimensional rectangular arrays, one for each
% component (e.g. {u,v,w}). User control allows representing vector data
% quiver arrows and color for magnitude, or as RGB color values. Quiver
% scaling is adjustable. If the field varies over time, provide
% four-dimensional rectangular arrays instead (with time varying along the
% last dimension). Optionally provide the spatial and temporal grid as a
% four-element cell array "xyzt", with the first element containing a
% vector of x locations, the second, y locations, the third, z locations,
% and the fourth, times (e.g. {x y z t}). Alternately, the four elements of
% "xyzt" can be four scalars indicating the grid sizes and time step (e.g.
% {dx dy dz dt}). The graphics handle of the resulting figure is provided
% as "hf". To save a movie, enter a file name ("movie name"), set the frame
% rate ("fps"), and click "play". Conversion from a grid made with ndgrid.m
% can be done via x_mesh = permute(x_nd,[2 1 3]). Inspired by vis4d.m by
% Joshua Stough, with some code re-used. 

% Written 30 January 2020 by Doug Kelley. 
% Update 16 February 2020: Can now play and save movies, accept vector
% field (or RGB images) as input, accept coordinates as inputs, allow
% toggling between vector and RGB plotting, allow changing the quiver
% scale, and provide figure handle as output. 

% Next: Make a similar visualization tool for 2D data. 

% -=- Constants -=-
figSize = [200 200]; % size of figure
minSpeedFactor = 2; % ignore drags slower than this
sliderPos = [0.645 0.045 0.2 0.025]; % position of time slider, normalized units
timeLabelPos = [0.48 0.048 0.115 0.025]; % position of time label, normalized units
timeTextPos = [0.595 0.05 0.05 0.025]; % position of time text box, normalized units
MovNameLabelPos = [0.5 0.02 0.095 0.025]; % position of movie name label, normalized units
MovNameTextPos = [0.595 0.02 0.205 0.025]; % position of movie name text box, normalized units
FrameRateLabelPos = [0.79 0.02 0.05 0.025]; % position of frame rate label, normalized units
FrameRateTextPos = [0.84 0.02 0.05 0.025]; % position of frame rate text box, normalized units
PlayPauseButtonPos = [0.9 0.02 0.05 0.025]; % position of play button, normalized units
VectorRgbGroupPos = [0.32 0.02 0.16 0.025]; % position of radio button group, normalized units
VectorButtonPos = [0 0 0.5 1]; % position of vector button, normalized units
RgbButtonPos = [0.5 0 0.5 1]; % position of rgb button, normalized units
QuiverScaleTextPos = [0.43 0.05 0.05 0.025]; % position of quiver scale text box, normalized
QuiverScaleLabelPos = [0.325 0.05 0.1 0.025]; % position of quiver scale label, normalized
xTextColor = [0.8 0.3 0.3]; % red-ish
yTextColor = [0.3 0.8 0.3]; % green-ish
zTextColor = [0.3 0.3 0.8]; % blue-ish
MaxQuiverCount = 30; % in each direction
qs=0.9; % initial scaling factor for quivers

% -=- Defaults -=-
xyzt_default = {1 1 1 1}; % grid size 1, time step 1
FrameRate_default = 10; % frames per second

% -=- Parse inputs and get set up -=-
if nargin<1
    error(['Usage: hf = ' mfilename '(im,[xyzt])'])
end
hasComps=false;
isRGB=false;
if iscell(im)
    if numel(im)==1
        im=im{1};
    elseif numel(im)==3
        imComps=im; % vector components
        hasComps=true;
        im=sqrt(im{1}.^2+im{2}.^2+im{3}.^2); % magnitude
    else
        error('Sorry, im must have either 1 or 3 components.')
    end
end
sz=size(im);
if numel(sz)==3
    sz(4)=1;
elseif numel(sz)>4 || numel(sz)<3
    error('Sorry, only 3D and 4D arrays are supported.')
end
if ~exist('xyzt','var') || isempty(xyzt)
    xyzt=xyzt_default;
end
if ~iscell(xyzt)
    error('Sorry, xyzt must be a cell array.')
end
if numel(xyzt)==3
    xyzt{4}=xyzt_default{4};
elseif numel(xyzt)~=4
    error('Size of xyzt does not match size of im.')
end
if all(cellfun(@numel,xyzt)==1)
    xyzt{1}=(0:sz(2)-1)*xyzt{1}; % coordinates, meshgrid-style
    xyzt{2}=(0:sz(1)-1)*xyzt{2};
    xyzt{3}=(0:sz(3)-1)*xyzt{3};
    xyzt{4}=(0:sz(4)-1)*xyzt{4};
elseif any(cellfun(@numel,xyzt)~=sz([2 1 3 4]))
    error('Size of xyzt does not match size of im.')
end
ind = isinf(im(:)) | isnan(im(:));
if any(ind)
    warning('Setting infinities and NaNs to zero.')
    im(ind)=0;
    if hasComps
        imComps{1}( isinf(imComps{1}) | isnan(imComps{1}) ) = 0;
        imComps{2}( isinf(imComps{2}) | isnan(imComps{2}) ) = 0;
        imComps{2}( isinf(imComps{2}) | isnan(imComps{2}) ) = 0;
    end
end
if hasComps
    ch1m=min(imComps{1}(:)); % in case the user wants RGB mode
    ch2m=min(imComps{2}(:));
    ch3m=min(imComps{3}(:));
    ch1r=range(imComps{1}(:));
    ch2r=range(imComps{2}(:));
    ch3r=range(imComps{3}(:));
    if ch1r==0  % avoid NaNs later
        ch1r=1; 
    end
    if ch2r==0  % avoid NaNs later
        ch2r=1; 
    end
    if ch3r==0  % avoid NaNs later
        ch3r=1; 
    end
end
minSpeed=minSpeedFactor*min( ...
    [mean(diff(xyzt{1})) mean(diff(xyzt{2})) mean(diff(xyzt{3}))]);
myxlim=[xyzt{1}(1) xyzt{1}(sz(2))]; % region extent, meshgrid-style. Calling "xlim" is slow.
myylim=[xyzt{2}(1) xyzt{2}(sz(1))];
myzlim=[xyzt{3}(1) xyzt{3}(sz(3))];
indx=round(sz(2)/2); % grid point nearest middle of region, meshgrid-style
indy=round(sz(1)/2);
indz=round(sz(3)/2);
indt=1; % first frame
if hasComps
    qxind=1:ceil(sz(2)/MaxQuiverCount):sz(2); % x indices of quivers
    qyind=1:ceil(sz(1)/MaxQuiverCount):sz(1); % y indices of quivers
    qzind=1:ceil(sz(3)/MaxQuiverCount):sz(3); % z indices of quivers
    maxlen=max(sqrt( ...
        reshape(imComps{1}(qyind,qxind,qzind,:),[],1).^2 ...
        + reshape(imComps{2}(qyind,qxind,qzind,:),[],1).^2 ...
        + reshape(imComps{3}(qyind,qxind,qzind,:),[],1).^2 )) ... % max magnitude
        / sqrt( ...
        mean(diff(xyzt{1})).^2 ...
        + mean(diff(xyzt{2})).^2 ...
        + mean(diff(xyzt{3})).^2 ); % diagonal grid size
end
lastPoint=[]; % initialize for later
t1=[]; % initialize for later

% -=- Plot everything for the first time -=-
hff=figure;
mypos=get(hff,'position');
mypos(3:4)=figSize;
set(hff,'position',mypos,'WindowButtonDownFcn',@buttonDownCallback, ...
    'color','w','WindowButtonUpFcn', @buttonUpCallback)

h=subplot(221);
hi=imagesc(xyzt{1},xyzt{2},squeeze(im(:,:,indz,indt)));
set(h(1),'xcolor',zTextColor,'ycolor',zTextColor,'nextplot','add')
if hasComps % plot quiver arrows for vector field
    [q1x,q1y]=meshgrid(xyzt{1}(qxind),xyzt{2}(qyind));
    hq=quiver(q1x(:),q1y(:), ...
        reshape(imComps{1}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
        reshape(imComps{2}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
        'k','autoscale','off');
end
hlx=plot(xyzt{1}(indx)*[1 1],myylim,'color',xTextColor);
hly=plot(myxlim,xyzt{2}(indy)*[1 1],'color',yTextColor);
xlabel('x')
ylabel('y')
ht=title(['z = ' num2str(xyzt{3}(indz),'%.3g') ...
    ' (drag up/down to change)'],'color',zTextColor);

h(2)=subplot(222);
hi(2)=imagesc(xyzt{3},xyzt{2},squeeze(im(:,indx,:,indt)));
set(h(2),'xcolor',xTextColor,'ycolor',xTextColor,'nextplot','add')
if hasComps % plot quiver arrows for vector field
    [q2z,q2y]=meshgrid(xyzt{3}(qzind),xyzt{2}(qyind));
    hq(2)=quiver(q2z(:),q2y(:), ...
        reshape(imComps{3}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
        reshape(imComps{2}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
        'k','autoscale','off');
end
hly(2)=plot(myzlim,xyzt{2}(indy)*[1 1],'color',yTextColor);
hlz=plot(xyzt{3}(indz)*[1 1],myylim,'color',zTextColor);
xlabel('z')
ylabel('y')
ht(2)=title(['x = ' num2str(xyzt{1}(indx),'%.3g') ...
    ' (drag up/down to change)'],'color',xTextColor);

h(3)=subplot(223);
hi(3)=imagesc(xyzt{1},xyzt{3},squeeze(im(indy,:,:,indt))');
set(h(3),'xcolor',yTextColor,'ycolor',yTextColor,'nextplot','add')
if hasComps % plot quiver arrows for vector field
    [q3x,q3z]=meshgrid(xyzt{1}(qxind),xyzt{3}(qzind));
    hq(3)=quiver(q3x(:),q3z(:), ...
        reshape(squeeze(imComps{1}(indy,qxind,qzind,indt))', ...
        [],1)/maxlen*qs, ...
        reshape(squeeze(imComps{3}(indy,qxind,qzind,indt))', ...
        [],1)/maxlen*qs, ...
        'k','autoscale','off');
end
hlx(2)=plot(xyzt{1}(indx)*[1 1],myzlim,'color',xTextColor);
hlz(2)=plot(myxlim,xyzt{3}(indz)*[1 1],'color',zTextColor);
xlabel('x')
ylabel('z')
ht(3)=title(['y = ' num2str(xyzt{2}(indy),'%.3g') ...
    ' (drag up/down to change)'],'color',yTextColor);

h(4)=subplot(224);
hs=slice(xyzt{1},xyzt{2},xyzt{3},im(:,:,:,indt), ...
    xyzt{1}(indx),xyzt{2}(indy),xyzt{3}(indz));
set(hs,'edgecolor','none');
axis tight
xlabel('x')
ylabel('y')
zlabel('z')
ht(4)=title(['t = ' num2str(xyzt{4}(indt),'%.3g')]);
hcb=colorbar;
hcb.Label.String='magnitude';
set(h(4),'nextplot','add','zdir','reverse')
if hasComps % plot quiver arrows for vector field
    hq(4)=quiver3(q1x(:),q1y(:),xyzt{3}(indz)*ones(size(q1x(:))), ...
        reshape(imComps{1}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
        reshape(imComps{2}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
        reshape(imComps{3}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
        'k','autoscale','off');
    hq(5)=quiver3(xyzt{1}(indx)*ones(size(q2y(:))),q2y(:),q2z(:), ...
        reshape(imComps{1}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
        reshape(imComps{2}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
        reshape(imComps{3}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
        'k','autoscale','off');
    hq(6)=quiver3(q3x(:),xyzt{2}(indy)*ones(size(q3x(:))),q3z(:), ...
        reshape(squeeze( ...
        imComps{1}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
        reshape(squeeze( ...
        imComps{2}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
        reshape(squeeze( ...
        imComps{3}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
        'k','autoscale','off');
end
hb=plot3( ...
    xyzt{1}(indx)*ones(1,5), ...
    [myylim(1) myylim(2) myylim(2) myylim(1) myylim(1)], ...
    [myzlim(2) myzlim(2) myzlim(1) myzlim(1) myzlim(2)], ...
    'color',xTextColor); % border lines
hb(2)=plot3( ...
    [myxlim(1) myxlim(2) myxlim(2) myxlim(1) myxlim(1)], ...
    xyzt{2}(indy)*ones(1,5), ...
    [myzlim(2) myzlim(2) myzlim(1) myzlim(1) myzlim(2)], ...
    'color',yTextColor);
hb(3)=plot3( ...
    [myxlim(1) myxlim(2) myxlim(2) myxlim(1) myxlim(1)], ...
    [myylim(2) myylim(2) myylim(1) myylim(1) myylim(2)], ...
    xyzt{3}(indz)*ones(1,5), ...
    'color',zTextColor);
myclim=[min(im(:)) max(im(:))];
set(h,'dataaspectratio',[1 1 1],'ydir','reverse','clim',myclim)
linkprop(h,'clim'); % one colorbar, so keep clim same for all axes
set([hlx hly hlz hb],'linewidth',2)

% -=- Make controls -=-
if sz(4)>1 % data includes more than 1 frame
    TimeLabel = uicontrol('Style','text','units','normalized', ...
        'Position',timeLabelPos,'FontSize',12,'FontWeight','bold', ...
        'String','t = ','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor','w');
    TimeText = uicontrol('Style','edit','units','normalized', ...
        'Position',timeTextPos,'String',num2str(xyzt{4}(1),'%.3g'), ...
        'Callback',@inputFromTimeText);
    TimeSlider = uicontrol('Style','slider','units','normalized', ...
        'Position',sliderPos,'Min',xyzt{4}(1),'Max',xyzt{4}(sz(4)), ...
        'Value',xyzt{4}(1),'SliderStep',[1/sz(4) 1/sz(4)+eps], ...
        'Callback',@inputFromTimeSlider);
    MovieNameLabel = uicontrol('Style','text','units','normalized', ...
        'Position',MovNameLabelPos,'FontSize',12,'FontWeight','bold', ...
        'String','movie name: ','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor','w');
    MovNameText = uicontrol('Style','edit','units','normalized', ...
        'Position',MovNameTextPos,'String','');
    FrameRateTextLabel = uicontrol('Style','text','units','normalized', ...
        'Position',FrameRateLabelPos,'FontSize',12,'FontWeight','bold', ...
        'String','fps: ','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor','w');
    FrameRateText = uicontrol('Style','edit','units','normalized', ...
        'Position',FrameRateTextPos,'String',FrameRate_default);
    PlayPauseButton = uicontrol('style','pushbutton', ...
        'units','normalized','position',PlayPauseButtonPos, ...
        'string','play','callback',@PlayPause);
else % if sz(4)>1
    set(ht(4),'visible','off')
end % if sz(4)>1
if hasComps
    VectorRgbGroup = uibuttongroup('units','normalized','position', ...
        VectorRgbGroupPos,'selectionchangedfcn',@VectorRgb);
    VectorButton = uicontrol(VectorRgbGroup,'style','radiobutton', ...
        'units','normalized','position',VectorButtonPos, ...
        'string','vector');
    RgbButton = uicontrol(VectorRgbGroup,'style','radiobutton', ...
        'units','normalized','position',RgbButtonPos, ...
        'string','rgb');
    QuiverScaleLabel = uicontrol('Style','text','units','normalized', ...
        'Position',QuiverScaleLabelPos,'FontSize',12,'FontWeight','bold', ...
        'String','quiver scale: ','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor','w');
    QuiverScaleText = uicontrol('Style','edit','units','normalized', ...
        'Position',QuiverScaleTextPos,'String',num2str(qs,'%.3g'), ...
        'Callback',@inputFromQuiverScaleText);
end % if hasComps
if nargout>0
    hf=hff;
end

% -=- Callback: When the user starts dragging -=-
    function buttonDownCallback(~,~)
        if gca==h(4)
            return % clicking on the 3D plot does nothing
        end
        p = get(gca,'CurrentPoint');
        lastPoint = [p(1); p(3)]; % get ready for dragging
        if lastPoint(1) < myxlim(1) || lastPoint(1) > myxlim(2) || ...
            lastPoint(2) < myylim(1) || lastPoint(2) > myylim(2)
            return % clicking outside the plotted region does nothing
        end
        set(hff,'WindowButtonMotionFcn',@dragCallback); % pay attention to dragging
    end % function buttonDownCallback(~,~)

% -=- Callback: When the user drags -=-
    function dragCallback(~,~)
        p = get(gca,'CurrentPoint');
        motionV = [p(1); p(3)] - lastPoint;
        if abs(motionV(2)) < minSpeed
            return
        end % Wait until the motion is noticeable.
        lastPoint = [p(1); p(3)];
        switch gca
            case h(1)
                [~,indz]=min(abs(xyzt{3}(indz)+motionV(2)-xyzt{3})); % drag down to increase z
                indz=min([max([indz 1]) sz(3)]); % don't allow coords outside the region
                if hasComps && isRGB
                    im1=cat(3, ...
                        (squeeze(imComps{1}(:,:,indz,indt))-ch1m)/ch1r, ... % rgb channels, with normalized intensity
                        (squeeze(imComps{2}(:,:,indz,indt))-ch2m)/ch2r, ...
                        (squeeze(imComps{3}(:,:,indz,indt))-ch3m)/ch3r );
                    im1(im1<0)=0;
                else
                    im1=squeeze(im(:,:,indz,indt));
                end
                set(hi(1),'cdata',im1) % update image
                if hasComps
                    set(hq(1), ... % update quivers
                        'udata',reshape( ...
                        imComps{1}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
                        'vdata',reshape( ...
                        imComps{2}(qyind,qxind,indz,indt),[],1)/maxlen*qs);
                end
                set(hlz(1),'xdata',xyzt{3}(indz)*[1 1]); % update lines on images
                set(hlz(2),'ydata',xyzt{3}(indz)*[1 1]);
                set(ht(1),'string',['z = ' num2str(xyzt{3}(indz),'%.3g') ...
                    ' (drag up/down to change)']) % update title
            case h(2)
                [~,indx]=min(abs(xyzt{1}(indx)-motionV(2)-xyzt{1})); % drag up to increase x
                indx=min([max([indx 1]) sz(2)]);
                if hasComps && isRGB
                    im2=cat(3, ...
                        (squeeze(imComps{1}(:,indx,:,indt))-ch1m)/ch1r, ...
                        (squeeze(imComps{2}(:,indx,:,indt))-ch2m)/ch2r, ...
                        (squeeze(imComps{3}(:,indx,:,indt))-ch3m)/ch3r );
                    im2(im2<0)=0;
                else
                    im2=squeeze(im(:,indx,:,indt));
                end
                set(hi(2),'cdata',im2)
                if hasComps
                    set(hq(2), ...
                        'udata',reshape( ...
                        imComps{3}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
                        'vdata',reshape( ...
                        imComps{2}(qyind,indx,qzind,indt),[],1)/maxlen*qs);
                end
                set(hlx,'xdata',xyzt{1}(indx)*[1 1]);
                set(ht(2),'string',['x = ' num2str(xyzt{1}(indx),'%.3g') ...
                    ' (drag up/down to change)'])
            case h(3)
                [~,indy]=min(abs(xyzt{2}(indy)+motionV(2)-xyzt{2})); % drag down to increase y
                indy=min([max([indy 1]) sz(1)]);
                if hasComps && isRGB
                    im3=cat(3, ...
                        (squeeze(imComps{1}(indy,:,:,indt))-ch1m)/ch1r, ...
                        (squeeze(imComps{2}(indy,:,:,indt))-ch2m)/ch2r, ...
                        (squeeze(imComps{3}(indy,:,:,indt))-ch3m)/ch3r );
                    im3(im3<0)=0;
                    im3=permute(im3,[2 1 3]);
                else
                    im3=squeeze(im(indy,:,:,indt))';
                end
                set(hi(3),'cdata',im3);
                if hasComps
                    set(hq(3), ...
                        'udata',reshape(squeeze( ...
                        imComps{1}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
                        'vdata',reshape(squeeze( ...
                        imComps{3}(indy,qxind,qzind,indt))',[],1)/maxlen*qs);
                end
                set(hly,'ydata',xyzt{2}(indy)*[1 1]);
                set(ht(3),'string',['y = ' num2str(xyzt{2}(indy),'%.3g') ...
                    ' (drag up/down to change)'])
            case h(4)
                return % dragging on the 3D plot does nothing
        end % switch gca
    end % function dragCallback(src,event)

% -=- Callback: When the user finishes dragging -=-
    function buttonUpCallback(~,~)
        set(hff, 'WindowButtonMotionFcn',''); % ignore mouse movements, for now
        set(hb(1),'xdata',xyzt{1}(indx)*ones(1,5)); % update border lines
        set(hb(2),'ydata',xyzt{2}(indy)*ones(1,5));
        set(hb(3),'zdata',xyzt{3}(indz)*ones(1,5));
        updateImagesAndSlices();
    end % function buttonUpCallback

% -=- Callback: When the user changes time using the slider bar -=-
    function inputFromTimeSlider(~,~)
        t1=get(TimeSlider,'value');
        changeTime();
    end % function inputFromTimeSlider(~,~)

% -=- Callback: When the user changes time using the text box -=-
    function inputFromTimeText(~,~)
        t1=str2double(get(TimeText,'string'));
        if isnan(t1) % user entered a non-number; replace with current time
            set(TimeText,'string',num2str(xyzt{4}(indt),'%.3g'));
            return
        end
        changeTime();
    end % function inputFromTimeText(~,~)

% -=- Callback: When the user changes quiver scale using the text box -=-
    function inputFromQuiverScaleText(~,~)
        qs1=str2double(get(QuiverScaleText,'string'));
        if isnan(qs1) % user entered a non-number; replace with current quiver scale
            set(QuiverScaleText,'string',num2str(qs,'%.3g'));
        else
            qs=qs1;
            updateImagesAndSlices();
        end
    end % function inputFromQuiverScaleText(~,~)

% -=- Callback: When user clicks PlayPauseButton -=-
    function PlayPause(~,~)
        if strcmpi(get(PlayPauseButton,'string'),'play')
            Play()
        else
            Pause()
        end
    end % function PlayPause(~,~)

% -=- Callback: When user chooses vector or rgb -=-
    function VectorRgb(~,event)
        if strcmpi(event.NewValue.String,'rgb')
            isRGB=true;
            set([hcb hq QuiverScaleText QuiverScaleLabel], ...
                'visible','off') % hide colorbar & quivers
        else
            isRGB=false;
            set([hcb hq QuiverScaleText QuiverScaleLabel], ...
                'visible','on') % show colorbar & quivers
        end
        updateImagesAndSlices();
    end % function VectorRgb(source,event)

% -=- Start playback -=-
    function Play(~)
        set(ht(1),'string',['z = ' num2str(xyzt{3}(indz),'%.3g')]) % update title
        set(ht(2),'string',['x = ' num2str(xyzt{1}(indx),'%.3g')])
        set(ht(3),'string',['y = ' num2str(xyzt{2}(indy),'%.3g')])
        FrameRate=str2double(get(FrameRateText,'string'));
        if isnan(FrameRate)
            FrameRate=FrameRate_default;
            set(FrameRateText,'string',num2str(FrameRate))
        end
        MovName=get(MovNameText,'string');
        if ~isempty(MovName) % save playback to disk
            indt0=indt; % save the time for later
            set([TimeLabel TimeText TimeSlider MovieNameLabel ...
                MovNameText FrameRateTextLabel FrameRateText ...
                PlayPauseButton],'visible','off')
            if hasComps
                set([VectorButton RgbButton VectorRgbGroup ...
                    QuiverScaleText QuiverScaleLabel],'visible','off')
            end
            [~,~,ext]=fileparts(MovName);
            if ~strcmpi(ext,'.avi')
                MovName=[MovName '.avi']; % append extension if necessary
            end
            vid=VideoWriter(MovName,'Motion JPEG AVI');
            vid.FrameRate=FrameRate;
            open(vid);
            for ii=1:sz(4)
                t1=xyzt{4}(ii);
                changeTime()
                pause(1/FrameRate);
                snap=getframe(hff);
                writeVideo(vid,snap.cdata);
            end % for ii=1:sz(4)
            close(vid);
            set([TimeLabel TimeText TimeSlider MovieNameLabel ...
                MovNameText FrameRateTextLabel FrameRateText ...
                PlayPauseButton],'visible','on')
            if hasComps
                set([VectorButton RgbButton VectorRgbGroup ...
                    QuiverScaleText QuiverScaleLabel],'visible','on')
            end
            t1=xyzt{4}(indt0);
            changeTime() % go back to the same time
            msgbox(['Saved ' MovName '.'])
        else % if ~isempty(MovName) - just play, don't save
            set(PlayPauseButton,'string','pause')
            set(hff,'WindowButtonDownFcn','','WindowButtonUpFcn','') % disable dragging
            set([TimeText MovNameText FrameRateText],'enable','off'); % disable controls
            if hasComps
                set([QuiverScaleText VectorButton RgbButton], ...
                    'enable','off'); % disable controls
            end
            set(TimeSlider,'callback',''); % disable slider while movie plays
            while strcmpi(get(PlayPauseButton,'string'),'pause')
                indt=indt+1;
                if indt>sz(4)
                    indt=1; % loop it
                end
                t1=xyzt{4}(indt);
                changeTime();
                pause(1/FrameRate);
            end % while strcmpi(get(PlayPauseButton,'string'),'pause')
        end % if ~isempty(MovName)
    end % function Play(~)

% -=- Pause playback -=-
    function Pause(~)
        set(PlayPauseButton,'string','play')
        set(hff,'WindowButtonDownFcn',@buttonDownCallback, ...
            'WindowButtonUpFcn', @buttonUpCallback); % re-enable dragging
        set([TimeText MovNameText FrameRateText],'enable','on'); % re-enable controls
        if hasComps
            set([QuiverScaleText VectorButton RgbButton],'enable','on'); % re-enable controls
        end
        set(TimeSlider,'callback',@inputFromTimeSlider); % re-enable slider
        set(ht(1),'string',['z = ' num2str(xyzt{3}(indz),'%.3g') ...
            ' (drag up/down to change)']) % reinstate reminders about dragging, in title
        set(ht(2),'string',['x = ' num2str(xyzt{1}(indx),'%.3g') ...
            ' (drag up/down to change)'])
        set(ht(3),'string',['y = ' num2str(xyzt{2}(indy),'%.3g') ...
            ' (drag up/down to change)'])
        t1=xyzt{4}(indt);
        changeTime() % just to clean up
        return
    end % function Pause(~)

% -=- Update controls and images and slices when time changes -=-
    function changeTime(~)
        [~,indt]=min(abs(t1-xyzt{4}));
        indt=min([max([indt 1]) sz(4)]); % don't allow times outside the range
        set(TimeText,'string',num2str(xyzt{4}(indt),'%.3g')) % update the text box
        set(TimeSlider,'Value',xyzt{4}(indt)); % update slider
        set(ht(4),'string',['t = ' num2str(xyzt{4}(indt),'%.3g')]);
        updateImagesAndSlices();
    end % function changeTime(~)

% -=- Update images, slices, and quivers -=-
    function updateImagesAndSlices(~)
        if hasComps && isRGB
            im1=cat(3, ...
                (squeeze(imComps{1}(:,:,indz,indt))-ch1m)/ch1r, ... % rgb channels, with normalized intensity
                (squeeze(imComps{2}(:,:,indz,indt))-ch2m)/ch2r, ...
                (squeeze(imComps{3}(:,:,indz,indt))-ch3m)/ch3r );
            im1(im1<0)=0;
            im2=cat(3, ...
                (squeeze(imComps{1}(:,indx,:,indt))-ch1m)/ch1r, ...
                (squeeze(imComps{2}(:,indx,:,indt))-ch2m)/ch2r, ...
                (squeeze(imComps{3}(:,indx,:,indt))-ch3m)/ch3r );
            im2(im2<0)=0;
            im3=cat(3, ...
                (squeeze(imComps{1}(indy,:,:,indt))-ch1m)/ch1r, ...
                (squeeze(imComps{2}(indy,:,:,indt))-ch2m)/ch2r, ...
                (squeeze(imComps{3}(indy,:,:,indt))-ch3m)/ch3r );
            im3(im3<0)=0;
        else
            im1=squeeze(im(:,:,indz,indt));
            im2=squeeze(im(:,indx,:,indt));
            im3=squeeze(im(indy,:,:,indt));
        end
        set(hi(1),'cdata',im1) % update images
        set(hi(2),'cdata',im2)
        if hasComps && isRGB
            set(hi(3),'cdata',permute(im3,[2 1 3]))
        else
            set(hi(3),'cdata',im3')
        end
        set(hs(1),'cdata',im2,'xdata',xyzt{1}(indx)*ones(sz([1 3]))); % update slices
        set(hs(2),'cdata',im3,'ydata',xyzt{2}(indy)*ones(sz([2 3])));
        set(hs(3),'cdata',im1,'zdata',xyzt{3}(indz)*ones(sz([1 2])));
        if hasComps
            set(hq(1), ... % update quivers
                'udata', ...
                reshape(imComps{1}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
                'vdata', ...
                reshape(imComps{2}(qyind,qxind,indz,indt),[],1)/maxlen*qs);
            set(hq(2), ...
                'udata', ...
                reshape(imComps{3}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
                'vdata', ...
                reshape(imComps{2}(qyind,indx,qzind,indt),[],1)/maxlen*qs);
            set(hq(3), ...
                'udata',reshape(squeeze( ...
                imComps{1}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
                'vdata',reshape(squeeze( ...
                imComps{3}(indy,qxind,qzind,indt))',[],1)/maxlen*qs);
            set(hq(4), ...
                'zdata',xyzt{3}(indz)*ones(size(q1x(:))), ...
                'udata', ...
                reshape(imComps{1}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
                'vdata', ...
                reshape(imComps{2}(qyind,qxind,indz,indt),[],1)/maxlen*qs, ...
                'wdata', ...
                reshape(imComps{3}(qyind,qxind,indz,indt),[],1)/maxlen*qs);
            set(hq(5), ...
                'xdata',xyzt{1}(indx)*ones(size(q2y(:))), ...
                'udata', ...
                reshape(imComps{1}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
                'vdata', ...
                reshape(imComps{2}(qyind,indx,qzind,indt),[],1)/maxlen*qs, ...
                'wdata', ...
                reshape(imComps{3}(qyind,indx,qzind,indt),[],1)/maxlen*qs);
            set(hq(6), ...
                'ydata',xyzt{2}(indy)*ones(size(q3x(:))), ...
                'udata',reshape(squeeze( ...
                imComps{1}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
                'vdata',reshape(squeeze( ...
                imComps{2}(indy,qxind,qzind,indt))',[],1)/maxlen*qs, ...
                'wdata',reshape(squeeze( ...
                imComps{3}(indy,qxind,qzind,indt))',[],1)/maxlen*qs);
        end % if hasComps
    end % function updateImagesAndSlices(~)

end % function slicei
