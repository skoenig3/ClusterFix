  function [behaviortime] = BehavioralIndex(behavind)
        %function turns indexes into times by parsing at breaks in continuity
        dind = diff(behavind);
        gaps =find(dind > 1);
        behaveind = zeros(length(gaps),50);
        if ~isempty(gaps)
            for gapind = 1:length(gaps)+1;
                if gapind == 1;
                    temp = behavind(1:gaps(gapind));
                elseif gapind == length(gaps)+1
                    temp = behavind(gaps(gapind-1)+1:end);
                else
                    temp = behavind(gaps(gapind-1)+1:gaps(gapind));
                end
                behaveind(gapind,1:length(temp)) = temp;
            end
        else
            behaveind =  behavind;
        end
        behaviortime = zeros(2,size(behaveind,1));
        for index=1:size(behaveind,1)
            rowfixind = behaveind(index,:);
            rowfixind(rowfixind == 0) = [];
            behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
        end
    end