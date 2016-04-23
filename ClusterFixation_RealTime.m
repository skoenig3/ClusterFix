% Real time Cluster Fix Emulator
tic
datafile = 'MP100712.1';
samprate = 5/1000; % in secs
variables = {'Dist','Vel','Accel','Rotation'};
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);

eyedat  = getEyeData(datafile,samprate);
fixationstats = cell(1,length(eyedat));

cndlop = 1; %template condition
%---Extract Paramters from Eye Traces---%
x = eyedat{cndlop}(1,:)*24+400; %converts data from [-400,400] to [0,800]
y = eyedat{cndlop}(2,:)*24+300; %converts data from [-300,300] to [0,600]
x = resample(x,samprate*1000,1);%up sample to 1000 Hz
y = resample(y,samprate*1000,1);%up sample to 1000 Hz
xss = filtfilt(flt,1,x);
yss = filtfilt(flt,1,y);
velx = diff(xss);
vely = diff(yss);
vel = sqrt(velx.^2+vely.^2);
accel = abs(diff(vel));
angle = 180*atan2(vely,velx)/pi;
vel = vel(200:end-200); %we do not care about the 1st or last 200 ms
accel = accel(200:end-199);
rot = zeros(1,length(xss)-2);
dist = zeros(1,length(xss)-2);
for a = 1:length(xss)-2;
    rot(a) = abs(angle(a)-angle(a+1));
    dist(a) = sqrt((xss(a)-xss(a+2)).^2 + (yss(a)-yss(a+2)).^2);
end
dist = dist(200:end-199);
rot(rot > 180) = rot(rot > 180)-180;
rot = 360-rot; %want rotation to be small so fixation values are all small
rot = rot(200:end-199);

points = [dist' vel' accel' rot'];
for ii = 1:size(points,2) %normalizes points to [0 1] by parameter
    thresh = mean(points(:,ii))+3*std(points(:,ii));%move outliers
    points((points(:,ii) > thresh),ii) = thresh;
    points(:,ii) = points(:,ii)-min(points(:,ii));
    points(:,ii) = points(:,ii)/max(points(:,ii));
end

%---Global Clustering---%
sil = zeros(1,5); %determines the number of clusters by comparing the ratio
%of intercluster and intracluster distances, faster mod of silhouette
for numclusts = 2:5
    T = kmeans(points(1:10:end,2:4),numclusts,'replicate',5);
    [silh] = InterVSIntraDist(points(1:10:end,2:4),T);
    sil(numclusts) = mean(silh);
end
sil(sil > 0.9*max(sil)) = 1;
numclusters = find(sil == max(sil));
T = kmeans(points,numclusters(end),'replicate',5);
meanvalues = zeros(max(T),size(points,2));
stdvalues = zeros(max(T),size(points,2));
for TIND = 1:max(T);
    tc = find(T == TIND);
    meanvalues(TIND,:) = mean(points(tc,:));
    stdvalues(TIND,:) = std(points(tc,:));
end

% determines fixation clusters by overlapping distributions in velocity
% and acceleration state space, here assumes gaussian distributions
[~, fixationcluster] = min(sum(meanvalues(:,2:3),2));
T(T == fixationcluster) = 100;
fixationcluster2 = find(meanvalues(:,2) < meanvalues(fixationcluster,2)...
    +3*stdvalues(fixationcluster,2));
fixationcluster2(fixationcluster2 == fixationcluster)= [];
for iii = 1:length(fixationcluster2);
    T(T == fixationcluster2(iii)) = 100;
end
T(T ~= 100) = 2;
T(T == 100) = 1;

fixationindexes =  find(T == 1)';
[fixationtimes] = BehavioralIndex(fixationindexes);
fixationtimes(:,(diff(fixationtimes,1) < 25)) = []; %25 ms duration threshold

% %---Local Re-Clusteirng---%
% notfixations = [];
% for ii = 1:length(fixationtimes);
%     %select points left and right of fixation for comparison
%     altind = fixationtimes(1,ii)-50:fixationtimes(2,ii)+50;
%     altind(altind < 1) = []; altind(altind > length(points)) = [];
%     POINTS = points(altind,:); %does not re-nomralize
%     sil = zeros(1,5);
%     for numclusts = 1:5
%         T = kmeans(POINTS(1:5:end,:),numclusts,'replicate',5);
%         [silh] = InterVSIntraDist(POINTS(1:5:end,:),T);
%         sil(numclusts) = mean(silh);
%     end
%     sil(sil > 0.9*max(sil)) = 1;
%     numclusters = find(sil == max(sil));  %it's dangerous to have too many clusters
%     T = kmeans(POINTS,ceil(median(numclusters)),'replicate',5);
%     rng = zeros(max(T),2*(size(POINTS,2)-1));
%     % determines fixation clusters by overlapping median values in velocity
%     % and acceleration state space, here we DO NOT assume gaussian distributions
%     % because there are not as many points and distributions rarely
%     % are normal
%     medianvalues = zeros(max(T),size(POINTS,2));
%     for TIND = 1:max(T);
%         tc = find(T == TIND);
%         if length(tc) == 1
%             rng(TIND,:) = ones(1,size(rng,2));
%             medianvalues(TIND,:) = POINTS(tc,:);
%         else
%             rng(TIND,:) = [max(POINTS(tc,1:end-1)) min(POINTS(tc,1:end-1))];
%             medianvalues(TIND,:) = median(POINTS(tc,:));
%         end
%     end
%     [~, fixationcluster] = min(sum(medianvalues(:,2:3),2));
%     T(T == fixationcluster) = 100;
%     fixationcluster2 = find((medianvalues(fixationcluster,2) < rng(:,2) & ...
%         (medianvalues(fixationcluster,2) > rng(:,5))) & (medianvalues(fixationcluster,3)...
%         < rng(:,3) & (medianvalues(fixationcluster,3) > rng(:,6))));
%     fixationcluster2( fixationcluster2 == fixationcluster)= [];
%     for iii = 1:length(fixationcluster2)
%         T(T == fixationcluster2(iii)) = 100;
%     end
%     T(T ~= 100) = 2;
%     T(T == 100) = 1;
%     notfixations = [notfixations altind(T == 2)];
% end

%%---Remove Points that are not fixations determing by Local Re-Clustering---%
% [~, ia, ~] = intersect(fixationindexes,notfixations);
% fixationindexes(ia) = [];
saccadeindexes = 1:size(points,1);
[~ , ~, ib] = intersect(fixationindexes,saccadeindexes);
saccadeindexes(ib) = [];

%---Consolidate & turn indexes into times---%
[saccadetimes] = BehavioralIndex(saccadeindexes);
[fixationtimes] = BehavioralIndex(fixationindexes);
tooshort = find(diff(fixationtimes,1) < 5); %removes accidental fixationtimes
notbehav = [];
for ii = 1:length(tooshort);
    notbehav = [notbehav fixationtimes(1,tooshort(ii)):fixationtimes(2,tooshort(ii))];
end
saccadeindexes = sort([saccadeindexes notbehav]);
tooshort = find(diff(saccadetimes,1) < 10); %10 ms duration threshold for saccades
notbehav = [];
for ii = 1:length(tooshort);
    notbehav = [notbehav saccadetimes(1,tooshort(ii)):saccadetimes(2,tooshort(ii))];
end
fixationindexes = sort([fixationindexes notbehav]);

[fixationtimes] = BehavioralIndex(fixationindexes);
fixationtimes(:,(diff(fixationtimes,1) < 25))= []; %25 ms duration threshold
fixationindexes = [];
for ii = 1:length(fixationtimes);
    fixationindexes = [fixationindexes fixationtimes(1,ii):fixationtimes(2,ii)];
end
[fixationtimes, fixations] = BehavioralIndexXY(fixationindexes,x,y);
saccadeindexes = 1:size(points,1);
[~ , ~, ib] = intersect(fixationindexes,saccadeindexes);
saccadeindexes(ib) = [];
toc

accel = accel + std(accel(fixationindexes)); 
accel = accel/.001^2;
vel = vel+ std(vel(fixationindexes));
vel = vel/.001;

scatter(vel(fixationindexes),accel(fixationindexes))

% statepsace = zeros(max(round(dist)),max(round(vel)),max(round(accel)),360);

vel = vel/5;
accel = accel/500;
va =round([vel(fixationindexes)' accel(fixationindexes)']);
statespace = zeros(1,max(round(va(:,2))));


mx = 0;
for i = max(va(:,1)):-1:1
    for ii = max(va(:,2)):-1:1
        temp = va(find(va(:,1) == i),:);
        if ~isempty(temp)
            mx = max(temp(:,2));
        end
        statespace(i) = mx;
    end
end
statespace = statespace(end:-1:1);
%%
statespace = [statespace(1)*ones(1,50) statespace statespace(end)*zeros(1,50)];
statespace = filtfilt(1/25*ones(1,25),1,statespace);
statespace = statespace(51:end-50);

statespace2D = sparse(round(max(vel)),round(max(accel)));
for i = 1:length(statespace);
    statespace2D(i,1:round(statespace(i))) = 1;
end
toc
%%
maxfixv = max(va(:,1))*2;
maxfixa = max(va(:,1))*2;
fixations = cell(1,length(eyedat));
extra = 0;
lag = 0;
for cndlop = 1:2:length(eyedat)
    type = NaN(1,length(x)); %1 fixation 0 sacccade
    time = 0;
    x = eyedat{cndlop}(1,:)*24+400;
    y = eyedat{cndlop}(2,:)*24+300;
    xpast = NaN(1,10);
    ypast = NaN(1,10);
    for tt = 1:length(x)-1;
        tic
        time = time + 0.005;
        xpast = [x(round(time/0.005)) xpast(1:end-1) ];
        ypast = [y(round(time/0.005)) ypast(1:end-1) ];
        if time > 0.05;
            velx = diff(xpast);
            vely = diff(ypast);
            vel = sqrt(velx.^2+vely.^2);
            accel = abs(diff(vel));
            accel = accel/.005^2/500;
            vel = vel/.005/5;
            velnow = round(mean(vel(end-7:end)));
            accelnow = round(mean(accel(end-7:end)));
            velnow(velnow > maxfixv) = maxfixv;
            accelnow(accelnow > maxfixa) = maxfixa;
            if full(statespace2D(velnow,accelnow)) == 1
                type(tt) = 1;
            else
                type(tt) = 0;
            end
        end
        tpass = toc;
        if tpass < 0.005
            extra = extra + 0.005-tpass;
        else
            lag = lag+0.005-tpass;
        end
    end
    
    type(isnan(type)) = [];
    saccind = find(type == 0);
    sacgap = find(diff(saccind) > 1);
    badsaccind = [];
    for i = 1:length(sacgap)
        if i == 1;
            temp = saccind(1:sacgap(i));
            if length(temp) < 5
                badsaccind = [badsaccind 1:sacgap(i)];
            end
        elseif i == length(sacgap);
            temp = saccind(sacgap(i)+1:length(saccind));
            if length(temp) < 5
                badsaccind = [badsaccind sacgap(i)+1:length(saccind)];
            end
        else
            temp = saccind(sacgap(i-1)+1:sacgap(i));
            if length(temp) < 5
                badsaccind = [badsaccind sacgap(i-1)+1:sacgap(i)];
            end
        end
    end
    saccind(badsaccind) = [];
    fixind = 10:length(x);
    [~, ia, ~] = intersect(fixind,saccind);
    fixind(ia) = [];
    fixgap = find(diff(fixind) > 1);
    badfixind = [];
    for i = 1:length(fixgap)
        if i == 1;
            temp = fixind(1:fixgap(i));
            if length(temp) < 10
                badfixind = [badfixind 1:fixgap(i)];
            end
        elseif i == length(fixgap);
            temp = fixind(fixgap(i)+1:length(fixind));
            if length(temp) < 10
                badfixind = [badfixind fixgap(i)+1:length(fixind)];
            end
        else
            temp = fixind(fixgap(i-1)+1:fixgap(i));
            if length(temp) < 10
                badfixind = [badfixind fixgap(i-1)+1:fixgap(i)];
            end
        end
    end
    fixind(badfixind) = [];
    length(find(diff(fixind) > 1))
    toc
    figure
    plot(x,y,'g');
    hold on
    fixgap = find(diff(fixind) > 1);
    for i = 1:length(fixgap)
        if i == 1;
            temp = fixind(1:fixgap(i));
            
        elseif i == length(fixgap);
            temp = fixind(fixgap(i)+1:length(fixind));
            
        else
            temp = fixind(fixgap(i-1)+1:fixgap(i));
        end
        plot(x(temp),y(temp),'r')
    end
    hold off
    pause(0.5)
    close
end