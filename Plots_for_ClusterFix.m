%% Useful Plotting Code
%% plots global points in 3D space-run to line 90
figure
plot3(points(:,2),points(:,3),points(:,4),'r.')
%% plot global points in 3D state space by global cluster-
% run to line 110 to get data by cluster but not yet classified
% run to line 126 for final fixation vs saccade clusters
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
%% Plot xy by global cluster
a = 1:length(x);
ia = 'rgbmcyk';
figure
hold on
plot(x(200:end-200),y(200:end-200),'k')
for TIND = 1:max(T);
    xss = a(T == TIND);
    yss = find(diff(xss) > 1);
    yss =  [1 yss];
    for ii = 2:length(yss);
        rng = xss(yss(ii-1)+1:yss(ii));
%         rng(rng > length(x)-200) = [];
%         rng(rng < 200) = [];
    rng = rng + 200;
        plot(x(rng),y(rng),ia(TIND))
    end
end

%% plot vel for local cluster
% run to line 157 if want by cluster
% run to line 189 if want fixation vs saccades
figure
hold on
xss = 'rgbmcyk';
vel = vel/max(vel);
plot(altind,vel(altind),'k')
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
% run to line 157 if want by cluster
% run to line 189 if want fixation vs saccades
xss = 'rgbmcyk';
figure
hold on
plot(x([altind(1)-100:altind(end)+80]+200),y([altind(1)-100:altind(end)+80]+200),'k')
for ii = 1:max(T);
    pp = altind(find(T == ii))+200;
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
%% plot local clusters in 3D space
% run to line 157 if want by cluster
% run to line 189 if want fixation vs saccades
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