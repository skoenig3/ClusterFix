%vel = 200/24*vel(1:2000);
%%
% screen_size = get(0, 'ScreenSize');
load('vel')
vidObj = VideoWriter('thresh.avi','Archival');
open(vidObj);
thresh = [0:40 mean(vel)+1/2*std(vel)];
figure
set(gca,'NextPlot','replaceChildren')
% set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
for i = 1:length(thresh);
    plot(vel,'r')
    ind = find(vel > thresh(i));
    dind = diff(ind);
    gaps = find(dind > 1);
    hold on
    if ~isempty(gaps)
        for gapind = 1:length(gaps)+1;
            if gapind == 1;
                plot(ind(1:gaps(gapind)),vel(ind(1:gaps(gapind))),'g');
            elseif gapind == length(gaps)+1
                plot(ind(gaps(gapind-1)+1:end),vel(ind(gaps(gapind-1)+1:end)),'g');
            else
                plot(ind(gaps(gapind-1)+1:gaps(gapind)),vel(ind(gaps(gapind-1)+1:gaps(gapind))),'g');
            end
        end
    else
        plot(vel,'g')
    end
    plot(1:length(vel),thresh(i),'k')
    hold off
    legend('Fixations','Saccades')
    xlabel('Time (ms)')
    ylabel('Velocity (dva/sec)')
    box off
    currFrame = getframe(gcf);
    if i == length(thresh)
        for ii = 1:66;
            writeVideo(vidObj,currFrame);
        end
    else
        for ii = 1:6;
            writeVideo(vidObj,currFrame);
        end
    end
end
close(vidObj)
disp('Finished Writing to Video')