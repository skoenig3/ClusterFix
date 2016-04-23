close all
%%
plot3(points(:,2),points(:,3),points(:,4),'r.')
%% plot points in 3D state space by global cluster
ia = 'rgbmcyk';
figure
hold on
for TIND = 1:max(T);
    plot3(points((T == TIND),2),points((T == TIND),3),points((T == TIND),4),...
        [ia(TIND)  '.'],'markersize',6);
end
xlabel('velocity')
ylabel('acceleration')
zlabel('rotation')
view(-7,32)
%% Plot xy by global cluster

x = eyedat{cndlop}(1,:)*24+400; %converts dva to pixel and data from [-400,400] to [0,800]
y = eyedat{cndlop}(2,:)*24+300; %converts dva to pixel and from [-300,300] to [0,600]
x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
x = resample(x,samprate*1000,1);%up sample to 1000 Hz
y = resample(y,samprate*1000,1);%up sample to 1000 Hz
xss = x(101:end-100); %remove buffer after filtering
yss = y(101:end-100); %remove buffer after filtering

a = 1:length(xss);
ia = 'rgbmcyk';
figure
hold on
plot(xss,yss,'k')
for TIND = 1:max(T);
    x = a(T == TIND);
    y = find(diff(x) > 1);
    y =  [1 y];
    for ii = 2:length(y);
        rng = x(y(ii-1)+1:y(ii));
        plot(xss(rng),yss(rng),ia(TIND))
    end
end
axis off
%% plot vel for local cluster
figure
hold on
xss = 'rgbmcyk';
vel = vel/(mean(vel)+3*std(vel));
plot(altind,vel(altind),'r')
for ii = 1:max(T);
    pp = altind(find(T == ii));
    ib = find(diff(pp) > 1);
    if isempty(ib);
        plot(pp,vel(pp),xss(ii));
    else
        for iii = 1:length(ib)+1;
            if iii == 1;
                plot(pp(1:ib(iii)),vel(pp(1:ib(iii))),[xss(ii)]);
            elseif iii == length(ib)+1;
                plot(pp(ib(end)+1:end),vel(pp(ib(end)+1:end)),xss(ii));
            else
                plot(pp(ib(iii-1)+1:ib(iii)),vel(pp(ib(iii-1)+1:ib(iii))),xss(ii));
            end
        end
    end
end

%% plot xy by local cluster
xss = 'rgbmcyk';
figure
hold on
plot(x([altind(1)-100:altind(end)+80]),y([altind(1)-100:altind(end)+80]),'k')
for ii = 1:max(T);
    pp = altind(find(T == ii));
    ib = find(diff(pp) > 1);
    if isempty(ib);
        plot(x(altind),y(altind),xss(ii));
    else
        for ia = 1:length(ib)+1;
            if ia == 1;
                plot(x(pp(1:ib(ia))),y(pp(1:ib(ia))),[xss(ii)]);
            elseif ia == length(ib)+1;
                plot(x(pp(ib(end)+1:end)),y(pp(ib(end)+1:end)),xss(ii));
            else
                plot(x(pp(ib(ia-1)+1:ib(ia))),y(pp(ib(ia-1)+1:ib(ia))),xss(ii));
            end
        end
    end
end
%%
ia = 'rgbmcyk';
figure
hold on
for iii = 1:max(T);
    pp = find(T == iii);
    plot3(POINTS(pp,2),POINTS(pp,3),POINTS(pp,4),[ia(iii)  '.']);
end
view(3)
xlabel('velocity')
ylabel('acceleration')
zlabel('rotation')
%%
close all
ia = 'rgbmcyk';
figure
hold on
plot(altind,vel(altind),'k')
for iii = 1:max(T);
    pp = find(T == iii);
    tc = find(diff(pp) > 1);
    if isempty(tc);
        plot(altind(pp),vel(altind(pp)),ia(iii));
    else
        for yss = 1:length(tc)+1;
            if yss == 1;
                plot(altind(pp(1:tc(yss))),vel(altind(pp(1:tc(yss)))),[ia(iii)]);
            elseif yss == length(tc)+1;
                plot(altind(pp(tc(end)+1:end)),vel(altind(pp(tc(end)+1:end))),ia(iii));
            else
                plot(altind(pp(tc(yss-1)+1:tc(yss))),vel(altind(pp(tc(yss-1)+1:tc(yss)))),ia(iii));
            end
        end
    end
end
xlabel('Time (ms)')
ylabel('Normalized Velocity')
%% MP100712.1 CND1 Fixation 17
figure
plot(x(altind(1)-500:altind(end)),y(altind(1)-500:altind(end)),'k')
hold on
plot(x(altind),y(altind),'r')