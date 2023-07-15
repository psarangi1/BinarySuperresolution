function amp_est=resolve_amp(amp1,thresh,set)
%%Automated Amplitude estimation
%amp1    : Initial set of potential amplitudes
%thresh  : match threshold
%set     : Combination technique - 0 - average, 1-local estimate, 2-best
%match from next sequence
%Output:
%amp_est : Estimated amplitude
%%
amp_est=amp1(:,1);
for i=2:size(amp1,2)
    if isempty(amp_est)
        break;
    end
    diff_mat=abs(amp1(:,i)-amp_est.');
    temp=[];
    for j=1:size(diff_mat,2)
         [match_amp,id_min]=min(diff_mat(:,j));
         %id_b,id_min
         if match_amp<=thresh
             if set==0
                temp=[temp;(amp_est(j)+amp1(id_min,i))/2];
             elseif set==1
                temp=[temp;amp_est(j)];
             else
                temp=[temp;amp1(id_min,i)];
             end
         end
    end
    amp_est=temp;
    %amp_est=uniquetol(amp_est,1e-2);
end
end

