function onlyspk_disp_oasis(spk_gt,spk_est,t_vec,th)
    id=find(spk_est>=th);
    hold on
    for i=1:size(spk_gt,1)
        line([spk_gt(i),spk_gt(i)],[0,0.2],'color','r')
        
    end
    for i=1:length(id)
        line([t_vec(id(i)),t_vec(id(i))],[0.3,spk_est(id(i))/max(spk_est)*0.2+0.5],'color','k')
    end
    count_gt=1;
    for i=2:length(spk_gt)
        if spk_gt(i)-spk_gt(i-1)<0.5
            count_gt=count_gt+1;
        else
            if count_gt>1
                text(spk_gt(i-1)+0.01,0.1,num2str(count_gt),'color','r','FontSize',10)
                count_gt=1;
            end
        end
    end
end