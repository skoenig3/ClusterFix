%% Plot disperion threshold, MST, and Cluster Fix
fltord = 60; %filter order
lowpasfrq = 30; %Low pass frequency cutoff
nyqfrq = 200 ./ 2; %nyquist frequency
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set008'};
tags = {'MP'};
algorithm = {'fixation','IDT-12','IDT-18','IDT-24','IDT-30','MST-NaN','MST-0.5','MST-0.75','MST-1.mat','MST-1.25'};

for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    statfiles = [];
    CND = 1;
    for alg = 1:length(algorithm);
        for mi = 1:length(matfiles.mat);
            str = strfind(matfiles.mat{mi},algorithm{alg});
            if ~isempty(str)
                for ii = 1:length(tags);
                    strt = strfind(matfiles.mat{mi},tags{ii});
                    if ~isempty(strt)
                        load(matfiles.mat{mi});
                        
                        figure
                        x = fixationstats{CND}.XY(1,:);
                        y = fixationstats{CND}.XY(2,:);
                        fixationtimes = fixationstats{CND}.fixationtimes;
                        
                        x = [x(20:-1:1) x x(end:-1:end-20)]; %add 20 ms buffer for filtering
                        y = [y(20:-1:1) y y(end:-1:end-20)]; %add 20 ms buffer for filtering
                        x = filtfilt(flt,1,x); %filter
                        y = filtfilt(flt,1,y); %filter
                        x = x(21:end-20); %remove buffer after filtering
                        y = y(21:end-20); %remove buffer after filtering
                        
                        vel = sqrt(diff(x).^2 + diff(y).^2);
                        vel = vel*200/24;%convert to dva  per sec since sampled at 200 Hz
                        
                        
                        subplot(1,2,1)
                        hold on
                        plot(5*(fixationtimes(1,1):fixationtimes(2,end)),vel,'g')
                        
                        fixationtimes(:,fixationtimes(2,:) < 1600) = [];
                        fixationtimes(:,fixationtimes(1,:) > 1790) = [];
                        
                        
                        for f = 1:size(fixationtimes,2)
                            plot(5*(fixationtimes(1,f):fixationtimes(2,f)),...
                                vel(fixationtimes(1,f):fixationtimes(2,f)),'r')
                        end
                        hold off
                        xlim([5*fixationtimes(1,1)+10 5*fixationtimes(2,end)-10])
                        ylabel('Velocity (dva/sec)')
                        xlabel('Time (ms)')
                        
                        subplot(1,2,2)
                        hold on
                        plot(x(1600:fixationtimes(2,end)),y(1600:fixationtimes(2,end)),'g')
                        
                        
                        for f = 1:size(fixationtimes,2)
                            plot(x(fixationtimes(1,f):fixationtimes(2,f)),...
                                y(fixationtimes(1,f):fixationtimes(2,f)),'r')
                        end
                        plot([450 474],[100 100],'k')
                        hold off
                        axis square
                        axis off
                        subtitle(algorithm{alg})
                    end
                end
            end
        end
    end
end
%%
cd('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 3 Cluster Fix vs Dispersion')

for fig = 1:length(algorithm)
    figure(fig)
    saveas(gcf,[algorithm{fig} '.fig'])
    saveas(gcf,[algorithm{fig} '.eps'])
end