function output = data_processing(serise_name, endtime, plot_tag, sampling_ratio)

rawdata = load(['rawdata/',serise_name]);
lambda_alpha = rawdata.lambda_alpha;
rawdata = rmfield(rawdata,'lambda_alpha');
sample_names = fieldnames(rawdata);


for i = 1:length(sample_names)
    temp = rawdata.(sample_names{i});
    t = temp(:,1)/1e9;
    t = t - t(1);
    exp_pl = temp(:,2);
    if sampling_ratio ~= 1
        sampling_index = 1:round(1/sampling_ratio):length(t);
        t = t(sampling_index);
        exp_pl = exp_pl(sampling_index);
    end
    offset = 0;%exp_pl(end);  % can change to zero
    end_index = length(t(t<endtime));
    alpha_abs = lambda_alpha(2,i);
    alpha_pl = lambda_alpha(4,i);
    lambda = lambda_alpha(1,i);
    if plot_tag == 1
        figure(999)
        hold on
        plot(t*1e9,exp_pl)
    end
    t = t(1:end_index);
    exp_pl = exp_pl(1:end_index);
    postdata(i) = v2struct(t,exp_pl,offset,alpha_abs,alpha_pl,lambda,end_index);
    if plot_tag == 1
        figure(999)
        hold on
        plot(t*1e9,exp_pl,'ro')
    end
end

if plot_tag == 1
    hold off
end

for i = 1:length(sample_names)
    postdata(i).name = sample_names{i};
    postdata(i).serise_name = serise_name;
end

output = postdata;

end