function onlyspk_disp(dt,spk_gt,spk_est)
    [d_mr2,dsec_mr2]=spk_distance(spk_gt,spk_est);
    hold on
    for i=1:size(dsec_mr2.okshift,2)
        line([spk_gt(dsec_mr2.okshift(1,i)),spk_gt(dsec_mr2.okshift(1,i))],[0,0.2],'color','r')
        line([spk_gt(dsec_mr2.okshift(1,i)),spk_est(dsec_mr2.okshift(2,i))],[0.3,0.5],'color','b')
        line([spk_est(dsec_mr2.okshift(2,i)),spk_est(dsec_mr2.okshift(2,i))],[0.6,0.8],'color','k')
    end
    for i=1:length(dsec_mr2.falsep)
        line([spk_est(dsec_mr2.falsep(i)),spk_est(dsec_mr2.falsep(i))],[0.6,0.8],'color','k')
        plot(spk_est(dsec_mr2.falsep(i)),0.7,'kx')
        
    end
    for i=1:length(dsec_mr2.miss)
        line([spk_gt(dsec_mr2.miss(i)),spk_gt(dsec_mr2.miss(i))],[0,0.2],'color','r')
        plot(spk_gt(dsec_mr2.miss(i)),0.1,'ro')
        
    end
    count_gt=1;
    for i=2:length(spk_gt)
        if spk_gt(i)-spk_gt(i-1)<0.25
            count_gt=count_gt+1;
        else
            if count_gt>1
                text(spk_gt(i-1)+0.01,0.1,num2str(count_gt),'color','r','FontSize',10)
                count_gt=1;
            end
        end
    end
    count_est=1;
    for i=2:length(spk_est)
        if spk_est(i)-spk_est(i-1)<0.25
            count_est=count_est+1;
        else
            if count_est>1
                text(spk_est(i-1)+0.01,0.7,num2str(count_est),'color','k','FontSize',10)
                count_est=1;
            end
        end
    end
end