screen_size = get(0, 'ScreenSize');
for i = 1:length(fixationstats);
    x = fixationstats{i}.XY;
    fixations = fixationstats{i}.fixations;
    fixationtimes = fixationstats{i}.fixationtimes;
    %     fixationtimes(:,fixationtimes(1,:) > 400) = [];
    figure
     set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    hold on
    plot(x(1,:),x(2,:),'g');
    for ii = 1:length(fixationtimes);
        plot(x(1,fixationtimes(1,ii):fixationtimes(2,ii)),...
            x(2,fixationtimes(1,ii):fixationtimes(2,ii)),'r')
        %pause(0.5)
%         plot(mean(x(1,fixationtimes(1,ii):fixationtimes(2,ii))),...
%             mean(x(2,fixationtimes(1,ii):fixationtimes(2,ii))),'*k')
%          plot(fixations(1,ii),fixations(2,ii),'b*')
    end
    %     legend('Saccades','Fixations')
    %     pause(2)
    %     close all
    box off
    axis off
    
end
%%
figure
hold on
plot(x,y,'g');
for ii = 1:length(fixationtimes);
    plot(x(fixationtimes(1,ii):fixationtimes(2,ii)),...
        y(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
    %     pause(0.5)
end
hold on
%%
figure
hold on
plot(200/24*vel,'g');%convert to from pixels to dva  per sec
for ii = 1:length(fixationtimes);
    plot(fixationtimes(1,ii):fixationtimes(2,ii),...
        200/24*vel(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
end
hold on
thresh =mean(200/24*vel)+std(200/24*vel);
plot(0:length(vel),thresh,'k')
xlabel('Time (ms)')
ylabel('Velocity (dva/sec)')
xlim([0 2100])
%%
figure
hold on
plot(1000*200/24*accel,'g');
for ii = 1:length(fixationtimes);
    plot(fixationtimes(1,ii):fixationtimes(2,ii),...
        1000*200/24*accel(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
end
hold on
thresh = 1000*200/24*(mean(accel)+std(accel));
plot(0:2100,thresh,'k')
xlabel('Time (ms)')
ylabel('Acceleration (dva/sec^2)')
xlim([0 2100])
ylim([0 7000])
%% Clustser Fix vs Thresholds in State Space
T = zeros(1,length(points));
T(fixationindexes)=1;
T(T==0) = 2;

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
view(2)

plot([0 1],[(mean(accel)+std(accel))/(mean(accel)+3*std(accel))...
    (mean(accel)+std(accel))/(mean(accel)+3*std(accel))],'k');
plot([(mean(vel)+std(vel))/(mean(vel)+3*std(vel))...
    (mean(vel)+std(vel))/(mean(vel)+3*std(vel))],[0 1],'k');
%%
figure
hold on
plot(dist/24,'g');
for ii = 1:length(fixationtimes);
    plot(fixationtimes(1,ii):fixationtimes(2,ii),...
        1/24*dist(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
end
hold on
thresh = (mean(dist)+std(dist))/24;
plot(0:2000,thresh,'k')
xlabel('Time (ms)')
ylabel('Distance (dva)')
xlim([0 2100])
ylim([0 1.2])
%%
figure
plot(x,y,'g');
hold on
for ii = 1:length(fixationtimes);
    plot(x(fixationtimes(1,ii):fixationtimes(2,ii)),y(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
end
%     xlim([90 300])
% ylim([90 450])
%%
cndlop = 1
x = fixationstats{cndlop}.XY(1,:);
y = fixationstats{cndlop}.XY(2,:);
fixationtimes = fixationstats{cndlop}.fixationtimes;

fltord = 8;
lowpasfrq = 30;
nyqfrq = 200 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);

velx = diff(x);
vely = diff(y);
vel = sqrt(velx.^2+vely.^2);
vel = filtfilt(flt,1,vel);
figure
plot(vel,'g');
hold on
for ii = 1:length(fixationtimes);
    plot(fixationtimes(1,ii):fixationtimes(2,ii),vel(fixationtimes(1,ii):fixationtimes(2,ii)),'r')
end
%% Plot fixations vs saccades for different algorithms
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};

for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    for CND = 1:2:72
        figure
        for mi = 1:length(matfiles.mat);
            str = strfind(matfiles.mat{mi},'-fixation'); %for MST-1 must put MST-1.mat or will grab MST-1.25 as well
            if ~isempty(str)
                for ii = 1:length(tags);
                    strt = strfind(matfiles.mat{mi},tags{ii});
                    if ~isempty(strt)
                        load(matfiles.mat{mi});
                        subplot(2,2,ii)
                        hold on
                        x = fixationstats{CND}.XY(1,:);
                        y = fixationstats{CND}.XY(2,:);
                        plot(x,y,'g')
                        fixationtimes = fixationstats{CND}.fixationtimes;
                        for f = 1:size(fixationtimes,2)
                            plot(x(fixationtimes(1,f):fixationtimes(2,f)),...
                                y(fixationtimes(1,f):fixationtimes(2,f)),'r')
                        end
                        axis off
                    end
                end
            end
        end
        pause
    end
end