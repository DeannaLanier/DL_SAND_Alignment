   function [locV,pks1,loc]=ReorderAlign_FindPeaks(X_intv,Mode,NumPeak,x4)
% find peaks inside confidence interval
sample=size(X_intv,1);
pks=cell(sample,1); 
loc=cell(sample,1);
pks1=cell(sample,1);

switch Mode
    case 'guide'
        for ind=1:sample
            [pks{ind},loc{ind}]=findpeaks(X_intv(ind,:),'SortStr','descend');
            if isempty(loc{ind})
                [pks{ind},loc{ind}]=max(X_intv(ind,:));
            end
            loc{ind}=sort(loc{ind}(1:NumPeak),2); %arrange horizontally
            pks1{ind}=X_intv(ind,loc{ind});
        end
    case 'rspd'
        for ind=1:sample
            [pks{ind},loc{ind}]=findpeaks(X_intv(ind,:),'SortStr','descend');
            if numel(loc{ind})<NumPeak
                loc{ind}(end+1:NumPeak)=ones(1,NumPeak-numel(loc{ind}));
            end 
            loc{ind}=sort(loc{ind}(1:NumPeak),2); %arrange horizontally
            pks1{ind}=X_intv(ind,loc{ind});
        end
    case 'multi'
        for ind=1:sample
        [pks{ind},loc{ind}]=findpeaks(X_intv{ind},'SortStr','descend');
            if numel(loc{ind})<NumPeak
                loc{ind}(numel(loc{ind})+1:NumPeak)=1*((x4(end)-x4(1))*(sample/2-ind)>0)+round(size(X_intv{ind},2))*((x4(end)-x4(1))*(sample/2-ind)<=0);%use the left/right boundary of X_intv as the missed peaks % 08/08/2021
            end
        loc{ind}=sort(loc{ind}(1:NumPeak),2);
        pks1{ind}=X_intv{ind}(loc{ind});
        end
end        
locV=cell2mat(loc);
pks1=cell2mat(pks1);
    end