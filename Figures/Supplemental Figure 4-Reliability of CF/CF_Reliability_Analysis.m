%% Analysis of Consistency Across several monkeys
addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');
dir = 'C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF';

cd(dir)

files = {'ReliabilityData_MP.mat','ReliabilityData_TT.mat','ReliabilityData_IW.mat'};
eyefiles = {'MP100712.1','TT100709.1','IW100709.1'};

alldata = cell(length(files),72);
acrossclassifications = cell(1,length(files));
for file = 1%:length(files);
    load(files{file})
    acrossclassifications{file} = allclassifications;
    eyetrace  = getEyeData(eyefiles{file},.005);
    for trial = 1:2:72
        classification = allclassifications{trial};
        if strcmpi(files{file},'ReliabilityData_MP.mat') && trial == 71;
            %for some reason this iteration didn't work right. Lots of errors with Matlab crashing and doing weird things
            classification(90,:) = [];
        elseif strcmpi(files{file},'ReliabilityData_TT.mat')
            %ran for 1000 iterations while others ran for 100 because 1000 kept crashing
            classification = classification(1:100,:);
        end
        potsacs = find(sum(classification) < size(classification,1));
        d = find(diff(potsacs) ~= 1);
       
        sactimes = [];
        sactimes = [potsacs(1);potsacs(d(1))];
        for dd = 1:length(d)-1;
            sactimes = [sactimes [potsacs(d(dd)+1);potsacs(d(dd+1))]];
        end
        sactimes = [sactimes [potsacs(d(end)+1);potsacs(end)]];
        
        sumclass = sum(classification);
        data = NaN(2,size(sactimes,2));
        for i = 1:size(sactimes,2)
            times = sactimes(1,i):sactimes(2,i);
            sac = sumclass(times);
            
            data(1,i) = 100*(size(classification,1)-min(sac))/size(classification,1); %reliability and using min cuz care about detection
            data(2,i) = sqrt((eyetrace{trial}(1,times(1))-eyetrace{trial}(1,times(end))).^2 + ...
                (eyetrace{trial}(2,times(1))-eyetrace{trial}(2,times(end))).^2); %saccade amplitude
        end
        alldata{file,trial} = data;
    end
end

figure
hold on
for file = 1:length(files);
    ad = cell2mat(alldata(file,:));
    plot(ad(2,:),ad(1,:),['.k'])
end
xlabel('Saccade Amplitude (dva)')
ylabel('% Consistently Labeled as a Saccade')
xlim([0 10])


bins = zeros(2,51);
bins(1,:) = 50:100;
for file = 1:length(files)
    for trial = 1:2:72;
        classification = acrossclassifications{file}{trial};
        if strcmpi(files{file},'ReliabilityData_MP.mat') && trial == 71;
            %for some reason this iteration didn't work right. Lots of errors with Matlab crashing and doing weird things
            classification(90,:) = [];
        elseif strcmpi(files{file},'ReliabilityData_TT.mat')
            %ran for 1000 iterations while others ran for 100 because 1000 kept crashing
            classification = classification(1:100,:);
        end
        classification = 100*sum(classification)/size(classification,1);
        for cl = 1:length(classification);
            temp = floor(classification(cl));
            if temp < 50
                temp = 100-temp;
            end
            binind = find(bins(1,:) == temp);
            bins(2,binind) = bins(2,binind)+1;
        end
    end
end

figure
percent_consistency = 100*bins(2,:)/sum(bins(2,:));
bar(bins(1,:),percent_consistency);
ylabel('% of Time points')
xlabel('% Consistency')
xlim([49 101])
box off
title('Percent of Points Conistent Across applications')


disp([num2str(percent_consistency(end)) '% of time points were classified the same 100% of the time'])
disp([num2str(sum(percent_consistency(50:end))) '% of time points were classified the same at least 99% of the time'])
disp([num2str(sum(percent_consistency(46:end))) '% of time points were classified the same at least 95% of the time'])
disp([num2str(sum(percent_consistency(41:end))) '% of time points were classified the same at least 90% of the time'])

%% Comparing how more replicates affect consistency results
addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');
dir = 'C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF';

cd(dir)

files = {'ReliabilityData_MP.mat','ReliabilityData_MP_10.mat','ReliabilityData_MP_25.mat'};
eyefiles = {'MP100712.1','MP100712.1','MP100712.1'};

for file = 1:length(files)
    alldata = cell(1,72);
    bins = zeros(2,51);
    bins(1,:) = 50:100;
    load(files{file})
    eyetrace = getEyeData(eyefiles{file},.005);
    for trial = 1:2:72;
        classification = allclassifications{trial};
        if trial == 71;
            %for some reason this iteration didn't work right on original
            %scan path so keep consistent across all
            classification(90,:) = [];
        end

        potsacs = find(sum(classification) < size(classification,1));
        d = find(diff(potsacs) ~= 1);
        
        sactimes = [];
        sactimes = [potsacs(1);potsacs(d(1))];
        for dd = 1:length(d)-1;
            sactimes = [sactimes [potsacs(d(dd)+1);potsacs(d(dd+1))]];
        end
        sactimes = [sactimes [potsacs(d(end)+1);potsacs(end)]];
        
        sumclass = sum(classification);
        data = NaN(2,size(sactimes,2));
        for i = 1:size(sactimes,2)
            times = sactimes(1,i):sactimes(2,i);
            sac = sumclass(times);
            
            data(1,i) = 100*(size(classification,1)-min(sac))/size(classification,1); %reliability and using min cuz care about detection
            data(2,i) = sqrt((eyetrace{trial}(1,times(1))-eyetrace{trial}(1,times(end))).^2 + ...
                (eyetrace{trial}(2,times(1))-eyetrace{trial}(2,times(end))).^2); %saccade amplitude
        end
        alldata{trial} = data;
        
        classification = 100*sum(classification)/size(classification,1);
        for cl = 1:length(classification);
            temp = floor(classification(cl));
            if temp < 50
                temp = 100-temp;
            end
            binind = find(bins(1,:) == temp);
            bins(2,binind) = bins(2,binind)+1;
        end
        
    end
    
    ad = cell2mat(alldata);
          
    figure
    plot(ad(2,:),ad(1,:),['.k'])
    xlabel('Saccade Amplitude (dva)')
    ylabel('% Consistently Labeled as a Saccade')
    xlim([0 10])
    title(files{file})
    
    figure
    percent_consistency = 100*bins(2,:)/sum(bins(2,:));
    bar(bins(1,:),percent_consistency);
    ylabel('% of Time points')
    xlabel('% Consistency')
    xlim([49 101])
    box off
    title(files{file})
    
    disp('')
    disp(files{file})
    disp([num2str(percent_consistency(end)) '% of time points were classified the same 100% of the time'])
    disp([num2str(sum(percent_consistency(50:end))) '% of time points were classified the same at least 99% of the time'])
    disp([num2str(sum(percent_consistency(46:end))) '% of time points were classified the same at least 95% of the time'])
    disp([num2str(sum(percent_consistency(41:end))) '% of time points were classified the same at least 90% of the time'])  
end
%% Results from k-means ++
addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');
dir = 'C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF';

cd(dir)

files = {'ReliabilityData_MP_PP.mat'};

for file = 1:length(files)
    bins = zeros(2,51);
    bins(1,:) = 50:100;
    load(files{file})
    for trial = 1:2:72;
        classification = allclassifications{trial};
        if trial == 71;
            %for some reason this iteration didn't work right on original
            %scan path so keep consistent across all
            classification(90,:) = [];
        end
        classification = 100*sum(classification)/size(classification,1);
        for cl = 1:length(classification);
            temp = floor(classification(cl));
            if temp < 50
                temp = 100-temp;
            end
            binind = find(bins(1,:) == temp);
            bins(2,binind) = bins(2,binind)+1;
        end
    end
    
    figure
    percent_consistency = 100*bins(2,:)/sum(bins(2,:));
    bar(bins(1,:),percent_consistency);
    ylabel('% of Time points')
    xlabel('% Consistency')
    xlim([49 101])
    box off
    title(files{file})
    
    disp('')
    disp(files{file})
    disp([num2str(percent_consistency(end)) '% of time points were classified the same 100% of the time'])
    disp([num2str(sum(percent_consistency(50:end))) '% of time points were classified the same at least 99% of the time'])
    disp([num2str(sum(percent_consistency(46:end))) '% of time points were classified the same at least 95% of the time'])
    disp([num2str(sum(percent_consistency(41:end))) '% of time points were classified the same at least 90% of the time'])
end
%% More thorough Analysis of MP data
addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');
dir = 'C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF';

cd(dir)

files = {'ReliabilityData_MP.mat'};
eyefiles = {'MP100712.1'};

fltord = 60; %filter order
lowpasfrq = 30; %Low pass frequency cutoff
nyqfrq = 200 ./ 2; %nyquist frequency
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter

alldata = cell(length(files),72);
acrossclassifications = cell(1,length(files));
for file = 1%:length(files);
    load(files{file})
    acrossclassifications{file} = allclassifications;
    eyetrace  = getEyeData(eyefiles{file},.005);
    for trial = 1:2:72
        classification = allclassifications{trial};
        if strcmpi(files{file},'ReliabilityData_MP.mat') && trial == 71;
            %for some reason this iteration didn't work right. Lots of errors with Matlab crashing and doing weird things
            classification(90,:) = [];
        elseif strcmpi(files{file},'ReliabilityData_TT.mat')
            %ran for 1000 iterations while others ran for 100 because 1000 kept crashing
            classification = classification(1:100,:);
        end
        potsacs = find(sum(classification) < size(classification,1));
        d = find(diff(potsacs) ~= 1);
        
                x = eyetrace{trial}(1,:);
                y = eyetrace{trial}(2,:);
                x = [x(20:-1:1) x x(end:-1:end-20)]; %add 20 ms buffer for filtering
                y = [y(20:-1:1) y y(end:-1:end-20)]; %add 20 ms buffer for filtering
                x = filtfilt(flt,1,x); %filter
                y = filtfilt(flt,1,y); %filter
                x = x(21:end-20); %remove buffer after filtering
                y = y(21:end-20); %remove buffer after filtering
        
        sactimes = [];
        sactimes = [potsacs(1);potsacs(d(1))];
        for dd = 1:length(d)-1;
            sactimes = [sactimes [potsacs(d(dd)+1);potsacs(d(dd+1))]];
        end
        sactimes = [sactimes [potsacs(d(end)+1);potsacs(end)]];
        
        sumclass = sum(classification);
        data = NaN(2,size(sactimes,2));
        for i = 1:size(sactimes,2)
            times = sactimes(1,i):sactimes(2,i);
            sac = sumclass(times);
            
            data(1,i) = 100*(size(classification,1)-min(sac))/size(classification,1); %reliability and using min cuz care about detection
            data(2,i) = sqrt((eyetrace{trial}(1,times(1))-eyetrace{trial}(1,times(end))).^2 + ...
                (eyetrace{trial}(2,times(1))-eyetrace{trial}(2,times(end))).^2); %saccade amplitude
            
            %commented code below allows you to plot examples of velocity
            %and position traces as classified. And displays % consistency
            
            if (data(2,i)< 2) && (data(2,i)> 1.75) %&& (data(1,i) > 97) %select size and % consistency
                if (length(x) >= (times(end)+21)) && ((times(1)-20) >= 1)
                    
                    velx = diff(x(times(1)-20:times(end)+21));
                    vely = diff(y(times(1)-20:times(end)+21));
                    vel  = sqrt(velx.^2+vely.^2);
                    vel = vel/(mean(vel)+3*std(vel));
                    figure
                    subplot(1,2,1)
                    hold on
                    plot(vel,'r')
                    plot(21:length(vel)-19,vel(21:end-19),'g')
                    plot([0 10],[0 0],'k')
                    hold off
                    axis off
                    subplot(1,2,2)
                    hold on
                    plot(x(times(1)-20:times(end)+20),y(times(1)-20:times(end)+20),'r')
                    plot(x(times),y(times),'g')
                    hold off
                    axis square
                    axis off
                    box off
                    subtitle([num2str(data(1,i)) '% ' num2str(data(2,i)) 'dva']);
                    close
                end
            end
        end
        alldata{file,trial} = data;
    end
end

figure
hold on
for file = 1:length(files);
    ad = cell2mat(alldata(file,:));
    plot(ad(2,:),ad(1,:),['.k'])
end
xlabel('Saccade Amplitude (dva)')
ylabel('% Consistently Labeled as a Saccade')
xlim([0 10])


bins = zeros(2,51);
bins(1,:) = 50:100;
for file = 1:length(files)
    for trial = 1:2:72;
        classification = acrossclassifications{file}{trial};
        if strcmpi(files{file},'ReliabilityData_MP.mat') && trial == 71;
            %for some reason this iteration didn't work right. Lots of errors with Matlab crashing and doing weird things
            classification(90,:) = [];
        elseif strcmpi(files{file},'ReliabilityData_TT.mat')
            %ran for 1000 iterations while others ran for 100 because 1000 kept crashing
            classification = classification(1:100,:);
        end
        classification = 100*sum(classification)/size(classification,1);
        for cl = 1:length(classification);
            temp = floor(classification(cl));
            if temp < 50
                temp = 100-temp;
            end
            binind = find(bins(1,:) == temp);
            bins(2,binind) = bins(2,binind)+1;
        end
    end
end

figure
percent_consistency = 100*bins(2,:)/sum(bins(2,:));
bar(bins(1,:),percent_consistency);
ylabel('% of Time points')
xlabel('% Consistency')
xlim([49 101])
box off
title('Percent of Points Conistent Across applications')

disp([num2str(percent_consistency(end)) '% of time points were classified the same 100% of the time'])
disp([num2str(sum(percent_consistency(50:end))) '% of time points were classified the same at least 99% of the time'])
disp([num2str(sum(percent_consistency(46:end))) '% of time points were classified the same at least 95% of the time'])
disp([num2str(sum(percent_consistency(41:end))) '% of time points were classified the same at least 90% of the time'])
%%
adbins = [0:0.5:9.5;0.5:0.5:10];
bincount = zeros(1,length(adbins));
ads = NaN(250,length(adbins));

for i = 1:size(ad,2);
    binind = find((adbins(1,:) < ad(2,i)) & (adbins(2,:) >= ad(2,i)));
    bincount(binind) = bincount(binind)+1;
    ads(bincount(binind),binind) = ad(1,i);
end

figure
hold on
x = 0:0.1:10;
xs = adbins(1,:);
sads = nanmean(ads);
params = sigm_fit(xs,sads,[],[],1);
ys = param(1)+(param(2)-param(1))./(1+exp(-param(4)*(param(3)-x)));

plot(xs,sads,'.k');
plot(x,ys);
%%
figure
bar(xs,nanmedian(ads))
box off
ylabel('Median % Consistency')
xlim([-0.5 10])
set(gca,'xTick',[])
%%
zads = zeros(size(ads));
zads(ads == 100) = 1;
szads = sum(zads);

figure
bar(xs,szads./bincount*100)
box off
ylabel('% of Saccades with 100% Consistency')
xlim([-0.5 10])
set(gca,'xTick',[])