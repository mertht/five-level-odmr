function plot_full_nv_odmr(rabi_freq, excitation_rate)
    % reads .MAT files from most recent call to "full_nv_odmr()"
    % and plots PL intensity curve.
    
    load('sweep.mat');
    load('photon_counts.mat');
    load('detune.mat');
    
    plot_title = ['ODMR with Rabi = ' num2str(rabi_freq, 3) ' and Excitation Rate = ' num2str(excitation_rate, 3)];


    figure
    plot(detune, photon_counts);
    xlabel('Microwave Excitation Frequency');
    ylabel('Photoluminescence');
    title(plot_title);
end