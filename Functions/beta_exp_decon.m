function x=beta_exp_decon(y_down,h,flip)
    %Greedy Beta expansion Decoding algorithm
    %Inputs: 
    %y_down: Downsampled measurement
    %h: FIR filter
    %flip: flag to flip the result 
    r=y_down;
    x=zeros(size(h)).';
    for i=1:length(h)
        if r/h(i)>=0.99
            x(i)=1;
            r=r-h(i);
        else
            x(i)=0;
        end
    end
    if flip==1
        x=flipud(x);
    end
end