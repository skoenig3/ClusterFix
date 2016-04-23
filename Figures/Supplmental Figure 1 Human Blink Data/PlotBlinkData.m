load('Pt011_blk1_eyedat.mat')
%%
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
%%
breaks = NaN(length(new_pos_time)-1,2);
for k = 1:length(new_pos_time)-1;
    breaks(k,1) =  new_pos_time{k}(end)+5;
    breaks(k,2) =  new_pos_time{k+1}(1)-5;
end
%%

figure
hold on
for k=6:10
    plot(new_pos_time{k},filtfilt(flt,1,new_pos_x{k}),'b');
    plot(new_pos_time{k},filtfilt(flt,1,new_pos_y{k}),'r');
    %     plot(new_pos_time{k},new_pos_x{k},'g');
    %     plot(new_pos_time{k},new_pos_y{k},'k');
end
legend('Horizontal eye trace','Vertical eye trace')
for k = 7:9
    fill([breaks(k,1) breaks(k,2) breaks(k,2)  breaks(k,1)],...
        [-400 -400 400 400],'k')
end
% xlim([new_pos_time{1}(1)-500 new_pos_time{5}(end)+500])
%%
vel = cell(1,5);
for k = 6:10
    velx = diff(filtfilt(flt,1,new_pos_x{k}));
    vely = diff(filtfilt(flt,1,new_pos_y{k}));
    vel{k} = sqrt(velx.^2 + vely.^2);
end

figure
hold on
for k=6:10
    plot(new_pos_time{k}(1:end-1),vel{k},'b')
end
xlabel('Time (seconds')
ylabel('Velocity')
for k = 7:9
    fill([breaks(k,1) breaks(k,2) breaks(k,2)  breaks(k,1)],...
        [0 0 15 15],'k')
end
%%
fixationstats = cell(1,5);
for k = 6:10
    fixs = ClusterFixMJ([new_pos_x{k}';new_pos_y{k}'],1/1000);
    fixationstats{k} = fixs;
end
%%
figure
hold on
for k = 6:10
    
    xss = filtfilt(flt,1,new_pos_x{k});
    yss = filtfilt(flt,1,new_pos_y{k})-600;
    
    plot(new_pos_time{k},xss,'g');
    plot(new_pos_time{k},yss,'g');
    
    fixtimes = fixationstats{k}.fixationtimes;
    start = new_pos_time{k}(1)-1;
    for f = 1:size(fixtimes,2);
        plot(start+fixtimes(1,f):start+fixtimes(2,f),xss(fixtimes(1,f):fixtimes(2,f)),'r')
        plot(start+fixtimes(1,f):start+fixtimes(2,f),yss(fixtimes(1,f):fixtimes(2,f)),'r')
    end
end
xlim([59000 62750])
plot([59000 59000],[0 -120],'k')
plot([breaks(8,1) breaks(8,2)],[200 200],'k')
plot([breaks(9,1) breaks(9,2)],[200 200],'k')
set(gca,'xTickLabel',{'0','500','1000','1500','2000','2500','3000','3500'})
%%
figure
hold on
for k = 6:10
    
    xss = new_pos_x{k};
    yss = new_pos_y{k};
    start = new_pos_time{k}(1)-1;
    
    oldfixtimes = fixationstats{k}.fixationtimes;
    oldsactimes = fixationstats{k}.saccadetimes;
    
    fixtimes = [];
    sactimes =[];
    for f = 1:size(oldfixtimes,2);
        if any((oldfixtimes(:,f) + start) < 59000) || any((oldfixtimes(:,f) + start) > 62750)
            if any((oldfixtimes(:,f) + start) < 59000)
                if ~all((oldfixtimes(:,f) + start) < 59000)
                    oldfixtimes((oldfixtimes(:,f) + start) < 59000) = 59000;
                end
            else
                if ~all((oldfixtimes(:,f) + start) > 62750)
                    oldfixtimes((oldfixtimes(:,f) + start) > 62750) = 62750;
                end
            end
        else
            fixtimes = [fixtimes oldfixtimes(:,f)];
        end
    end
    for f = 1:size(oldsactimes,2);
        if any((oldsactimes(:,f) + start) < 59000) || any((oldsactimes(:,f) + start) > 62750)
            if any((oldsactimes(:,f) + start) < 59000)
                if ~all((oldsactimes(:,f) + start) < 59000)
                    oldsactimes((oldsactimes(:,f) + start) < 59000) = 59000;
                end
            else
                if ~all((oldsactimes(:,f) + start) > 62750)
                    oldsactimes((oldsactimes(:,f) + start) > 62750) = 62750;
                end
            end
        else
            sactimes = [sactimes oldsactimes(:,f)];
        end
    end
    
    for f = 1:size(sactimes,2);
        plot(xss(sactimes(1,f):sactimes(2,f)),yss(sactimes(1,f):sactimes(2,f)),'g')
    end
    for f = 1:size(fixtimes,2);
        plot(xss(fixtimes(1,f):fixtimes(2,f)),yss(fixtimes(1,f):fixtimes(2,f)),'r')
    end
end
hold off
%%
figure
hold on
for k = 6:10
    
    plot(new_pos_time{k}(1:end-1),vel{k},'g');
    
    fixtimes = fixationstats{k}.fixationtimes;
    start = new_pos_time{k}(1);
    for f = 1:size(fixtimes,2);
        fixt = fixtimes(1,f):fixtimes(2,f);
        fixt(fixt > length(vel{k})) = length(vel{k});
        plot(fixt+start-1,vel{k}(fixtimes(1,f):fixtimes(2,f)),'r')
    end
end
xlabel('Time (seconds')
ylabel('Normalized Velocity')
xlim([59000 62750])
plot([breaks(8,1) breaks(8,2)],[14 14],'k')
plot([breaks(9,1) breaks(9,2)],[14 14],'k')
set(gca,'xTickLabel',{'0','500','1000','1500','2000','2500','3000','3500'})
