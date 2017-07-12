function plot_cw_sensitivity()
    % Generates contour plot of magnetic field sensitivity vs
    % optical pumping and rabi frequency.
    
    % simulation parameters
    samples = 20; % number of data points
    opt_pow = linspace(0.01*10^9, 10^9, samples)'; % values of optial power
    rabi_freqs = linspace(100*10^3, 10 * 10^6, samples)'; % values of rabi frequencies
    
    % declare arrays for data to be plotted
    contrast_plot = zeros(samples, samples);
    photon_count_plot = zeros(samples, samples);
    sensitivity_plot = zeros(samples, samples);
    line_width_plot = zeros(samples, samples);

    % create plots
    for a = 1:length(opt_pow)
        for b = 1:length(rabi_freqs)
            [FWHM, contrast, photon_count] = single_nv_odmr(1 / (1 * 10^(-6)), 3 * 10^9, opt_pow(a), rabi_freqs(b), 0);
            contrast_plot(a, b) = contrast;
            photon_count_plot(a, b) = photon_count; % calculate total photoluminescence (per pixel)
            line_width_plot(a, b) = FWHM;
            deltaB = FWHM / (contrast * sqrt(photon_count)); % measure of magnetic sensitivity
            sensitivity_plot(a, b) = deltaB;
        end
    end
    
    % save data to csv files
    csvwrite('line_width_plot.csv', line_width_plot);
    csvwrite('contrast_plot.csv', contrast_plot);
    csvwrite('photon_count_plot.csv', photon_count_plot);
    csvwrite('sensitivity_plot.csv', sensitivity_plot);
    
    % plot sensitivity
    figure
    surf(opt_pow, rabi_freqs, sensitivity_plot);
    set(gca, 'zscale', 'log')
    title('Log Magnetic Field Sensitivity vs. Optical Power and Rabi Frequency')
    xlabel('Optical Pump Rate')
    ylabel('Rabi Frequency')
    
    % plot line width
    figure
    surf(opt_pow, rabi_freqs, line_width_plot);
    title('Line Width vs. Optical Power and Rabi Frequency')
    xlabel('Optical Pump Rate')
    ylabel('Rabi Frequency')

    % plot contrast
    figure
    surf(opt_pow, rabi_freqs, contrast_plot);
    %set(gca, 'zscale', 'log')
    title('Contrast vs. Optical Power and Rabi Frequency')
    xlabel('Optical Pump Rate')
    ylabel('Rabi Frequency')

    % plot photon count
    figure
    surf(opt_pow, rabi_freqs, photon_count_plot);
    %set(gca, 'zscale', 'log')
    title('Photon Count vs. Optical Power and Rabi Frequency')
    xlabel('Optical Pump Rate')
    ylabel('Rabi Frequency')
    
end