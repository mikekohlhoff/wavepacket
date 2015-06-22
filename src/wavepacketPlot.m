%%
% Plot wavefunction probability density
clear; clc; clf
%vel = 0.00014;   % initial velocity
%nx = 100;        % angular points
%ny = 45;        % radial points
%dt = 1;         % time step
%plotstep = 3000; % plot steps
%d0 = 30;       % initial surface distance
%aba = d0;       % absorbing boundary
%grid = d0+30;   % grid boundary

% Threshold for plotting wavefunction {0,1}
limc = 0.9;

%%
% Load in the wavefunction
[FileName,PathName] = uigetfile('wavefunction.*','Select the wavefunction file','*');
wf = load([PathName FileName]);
params = wf(1,:);
wf = wf(2:end,:);

%%
ny = params(1);
nx = params(2);
tstep = params(3);
aba = params(4);

%%
n = nx*ny;

% wavefunction
nt = length(wf)/n;
x = wf(1:n,1); y = wf(1:n,2); dd = wf(:,4); wf = wf(:,3);
wmax = max(wf); wmin = min(wf);
x = reshape(x,nx,ny); y = reshape(y,nx,ny);
dd = reshape(dd,n,[]); dd = dd(1,:)';

% axis for plotting
xmax = max(x(:));
xmax = xmax*cos(pi/4); xmin = -xmax;
rmax = ceil(sqrt(max(max(x.^2+y.^2))));

% absorbing boundary 'ring'
abx  =[0:aba]';
rabs = aba^2; aby = sqrt(rabs-abx.^2);
abx = [flipud(abx); -abx; flipud(-abx); abx];
aby = [flipud(aby); aby; flipud(-aby); -aby];

% grid ring
gridx = [0:rmax]';
rgrids = rmax^2; gridy = sqrt(rgrids-gridx.^2);
gridx = [flipud(gridx); -gridx; flipud(-gridx); gridx];
gridy = [flipud(gridy); gridy; flipud(-gridy); -gridy];

%%
%2D plot
clf;
w2 = wf/min(wf);  % set minimum to 1 so can take log10
w2 = log10(w2);
w2 = w2./max(w2); % normalise to max

iw2 = find(w2<limc); w2(iw2)=limc; % filter off wavefunction thats too small
wmax2 = max(w2); wmin2 = min(w2);
w2 = reshape(w2,nx,ny,[]);

h2 = pcolor([x;flipud(-x);x(1,:)],[y;flipud(y);y(1,:)],...
    [w2(:,:,1);flipud(w2(:,:,1));w2(1,:,1)]);hold on;shading interp; %plot wf

%plot(abx,aby,'y:','linewidth',1.5); % plot absorbing boundary

%plot(gridx,gridy,'b','linewidth',1.5); %plot grid edge

wall = patch([-rmax,-rmax,rmax,rmax],[-rmax,-dd(1),-dd(1),-rmax],[1 1 1],...
    'facealpha',0.5,'edgecolor','none'); % plot surface

txt = text(-rmax/5,rmax*0.8,['distance=',num2str(dd(1)),'a_{0}'],...
    'color','w','fontweight','bold'); % display current distance

core = plot(0,0,'o','MarkerFaceColor','b','MarkerSize',8,'markeredgecolor','none');

axis off; ([-rmax rmax -rmax rmax]); %axis square;
colormap('hot'); caxis([wmin2,wmax2]);

set(gca,'color',[0 0 0]);
set(gcf,'color',[0 0 0]);

% setup figure position and control buttons for video
set(gcf,'Units','normalized')
%set(gcf,'position',[0.1 0.2 1 0.5]);
h = uicontrol('Style', 'pushbutton', 'String', 'play',...
    'Position', [20 200 50 35], 'Callback', 'uiresume(gcbf)');

h = uicontrol('Style', 'pushbutton', 'String', 'pause',...
    'Position', [20 150 50 35], 'Callback', 'uiwait(gcf)');

h = uicontrol('Style', 'pushbutton', 'String', 'replay',...
    'Position', [20 100 50 35], 'Callback', 'uiresume(gcbf);iback =1;');

h = uicontrol('Style', 'pushbutton', 'String', 'stop',...
    'Position', [20 50 50 35], 'Callback', 'uiresume(gcbf);ibreak =1;');

%mov = avifile('WPP_video.avi')

uiwait(gcf);
iback = 0; ibreak = 0; it = 0;
while iback == 0
    it = it +1;
    if it > nt
        uiwait(gcf);
    else
        dist = dd(it);
        
        set(txt,'string',['distance = ',num2str(dist),'a_{0}']);
        
        set(h2,'cdata',([w2(:,:,it);flipud(w2(:,:,it));w2(1,:,it)]));
        set(wall,'ydata',[-rmax,-dist,-dist,-rmax])
        
        pause(0.05)
        drawnow
    end
    
    if (iback == 1) % set replay option
        iback = 0;
        it = 1;
    end
    
    if (ibreak == 1) % set stop video option
        break
    end
    
    %    if rem(it,4)==0 % set output steps to record movie
    %        F = getframe(gcf);
    %        mov = addframe(mov,F)
    %    end
    
end
%mov = close(mov)

%% 3d Plot

% view
ang = 95; el = 40; ax = [xmin xmax xmin xmax 0 1];

w1 = wf./wmax;
w1 = reshape(w1,nx,ny,[]);

clf;

h(1) = axes('Position',[0 0 1 1]); hold on;
xm = xmax+100;
surf([-xm -xm; xm xm],[-xm xm;-xm xm],[0,0;0,0],'facecolor',[0  0 0]);
axis(ax); axis off; view(ang,el)

h(2) = axes('Position',[0 0 1 1]); hold on;
pp = surf([x;flipud(-x);x(1,:)],[y;flipud(y);y(1,:)],...
    [w1(:,:,1);flipud(w1(:,:,1));w1(1,:,1)],'edgecolor','k',...
    'edgealpha',0.15,'facecolor','interp'); hidden off;

% surface
x1 = [xmin xmin xmax xmax]; y1 = -dd(1)*ones(1,4); z1 = [0 1 1 0];
wall = fill3(x1,y1,z1,'g','facealpha',0.5,'linestyle','--','edgecolor','g');

core = plot3(0,0,.01,'o','MarkerFaceColor','b','MarkerSize',8,'markeredgecolor','none');

txt = text(xmin,0,0.8,['Distance = ',num2str(dd(1)),'a_{0}'],...
    'color','w','fontweight','bold','fontsize',11);

set(gca, 'color', 'w'); set(gcf,'color',[0.5 0.5 0.5]);
axis(ax); axis off; view(ang,el)
xlabel('p');ylabel('z');colormap('hot');


set(gcf,'Units','normalized')
%set(gcf,'position',[0.1 0.2 1 0.5]);
h = uicontrol('Style', 'pushbutton', 'String', 'play',...
    'Position', [500 200 50 35], 'Callback', 'uiresume(gcbf)');

h = uicontrol('Style', 'pushbutton', 'String', 'pause',...
    'Position', [500 150 50 35], 'Callback', 'uiwait(gcf)');

h = uicontrol('Style', 'pushbutton', 'String', 'replay',...
    'Position', [500 100 50 35], 'Callback', 'uiresume(gcbf);iback =1;');

h = uicontrol('Style', 'pushbutton', 'String', 'stop',...
    'Position', [500 50 50 35], 'Callback', 'uiresume(gcbf);ibreak =1;');

%mov = avifile('WPP_video3d.avi')

uiwait(gcf);
iback = 0; ibreak = 0; it = 0;
while iback == 0
    it = it +1;
    if it > nt
        uiwait(gcf);
    else
        dist = dd(it);
        
        set(txt,'string',['distance = ',num2str(dist),'a_{0}']);
        
        set(pp,'zdata',([w1(:,:,it);flipud(w1(:,:,it));w1(1,:,it)]));
        set(wall,'ydata',-dd(it)*ones(1,4))
        
        pause(0.05)
        drawnow
    end
    
    if (iback == 1) % set replay option
        iback = 0;
        it = 1;
    end
    
    if (ibreak == 1) % set stop video option
        break
    end
    
    %    if rem(it,4)==0 % set output steps to record movie
    %        F = getframe(gcf);
    %        mov = addframe(mov,F)
    %    end
end
%mov = close(mov)
