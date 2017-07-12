function [FWHM, contrast, photon_count] = full_nv_odmr()
    % Simulates ODMR on an NV-center. Use for qualitatively determining the
    % effect of optical power, microwave power, and dephasing rates on
    % magnetic field sensitivity.
    
    clear all
    close all
    clc


    %% EMPIRICAL CONSTANTS
    
    nv_density          = 1 * 10^17 / (10^(-2))^3;      % average number of nv per m^3
    sample_width        = 10 * 10^(-6);                 % side length of diamond sample (m)
    sample_thickness    = 100 * 10^(-9);                % depth dimension of NV diamond sample (m)
    collected_photons   = 0.01;                         % experimental ratio of collected photons
    beam_sigma          = 3 * 10^(-6)/(2.35);           % std dev of gaussian beam (m)
    D0                  = 2.87 * 10^9;                  % zero-field splitting (Hz)
    G2star              = 1 / (1 * 10^(-6));            % spin-spin dephasing rate (1/s)
    rabi_frequency      = 2 * pi * 10^5;                % rabi frequency
    pixel_dim           = 1 * 10^(-7);                  % length dimension of sensor pixel (after binning) (m)
    excitation_rate     = 10^7;                         % pump rate of laser (Hz)
    
    
    %% SIMULATION PARAMETERS (change these for perfomance purposes)
    
    max_sim_res         = 75;                           % image cannot exceed this resolution length
    max_detune          = .03 * 10^9;                   % maximum deviation from w1
    freq_samples        = 50;                          % number of MW frequency slices to test
    plot_images         = 1;                            % "1" to plot expected captured image
    plot_intensity      = 1;                            % "1" to plot a PL profile over MW freq. sweep
    
    
    %% CHECK PARAMETERS
    
    res = min(sample_width / pixel_dim, max_sim_res);   % number of pixels along one side length
    min_pixel_size = sample_width / res; % now we have the size dim. of a simulation pixel
    
    % get scaling factor from a single NV to a pixel
    intensity_factor = (min_pixel_size^2) * sample_thickness * nv_density * collected_photons;

    if res == max_sim_res % warn user
        disp('defaulting to max_sim_res');
    end
    
    if max_detune >= D0
        error('max detune is too large')
    end
    
    
    %% CREATE NV DIAMOND IMAGE AND GAUSSIAN BEAM (X and Y coordinates)
    
    X = linspace(-sample_width/2, sample_width/2, res);
    Y = X;
    image = zeros(length(X), length(Y));
    
    Z = radial_gaussian(X, Y, beam_sigma);
    peak = max(max(Z));
    beam = (Z / peak) * excitation_rate;
    view(2);
        
    
    %% GENERATE BROADENING DISTRIBUTION (inhomogenous resonant frequency)
    
    % i.e. due to miniscule local changes in magnetic field
    gamma(1:res, 1:res) = 25 / 10^(-3); % 1/T2
    w1_sigma = 1 / (pi * G2star) / 2.35;
    w1 = normrnd(D0, w1_sigma, res, res); % sample splitting from Gaussian
    
    detune = linspace(D0 - max_detune, D0 + max_detune, freq_samples);
    sweep = zeros(freq_samples, res, res); % array of images (one for each tested MW frequency)
    photon_counts = zeros(1, length(detune)); % averaged PL value from center of image
    

    %% GET BACKGROUND PHOTON COUNT (MW off)
    
    mask = zeros(res, res);
    for a = 1:length(X)
        for b = 1:length(Y)
            % calculate and scale result of single nv odmr
            PL = single_nv_odmr(gamma(a, b), w1(a, b), 0, beam(a, b), D0);
            mask(a, b) = PL * intensity_factor;
        end
    end
    s2 = surf(X, Y, mask);
    set(s2,'LineStyle','none')
    title('Mask (no microwave excitation)')
    view(2);

    
    %% PERFORM FREQUENCY SWEEP
    
    for index = 1:length(detune)
        % perform cw ODMR on NV centers (pixel by pixel)
        for a = 1:length(X)
            for b = 1:length(Y)
                % calculate and scale result of single nv odmr
                PL = single_nv_odmr(gamma(a, b), w1(a, b), rabi_frequency, beam(a, b), detune(index));
                image(a, b) = PL * intensity_factor;
            end
        end
        message = strcat('finished image  ', int2str(index));
        disp(message)
        net = image - mask; % get difference images
        sweep(index, :, :) = net;
        
        photon_counts(index) = get_pl(net); % get PL from the central region of the image
        if plot_images
            figure
            s2 = surf(X, Y, net);
            set(s2,'LineStyle','none')
            view(2);
            colorbar
        end
    end
    
    
    %% PLOT RESULTS
    
    save('sweep.mat', 'sweep');
    save('photon_counts.mat', 'photon_counts');
    save('detune.mat', 'detune');

    if plot_intensity
        plot_full_nv_odmr(rabi_freqency, excitation_rate);
    end
    
    
    %% FIND DISTRIBUTION PARAMETERS
    
    photon_count = max(photon_counts);
    norm_pl = photon_counts / photon_count; % normalize PL
    
    contrast = 1 - min(norm_pl);
    peak_left = -1; % index of left FWHM
    peak_right = -1; % index of right FWHM
    found_left = 0; % boolean flag
    
    % scan for peak width
    for n = 1:length(norm_pl)
        if (norm_pl(n) < 1 - contrast/2) && (~found_left)
            peak_left = n;
            found_left = 1;
        end
        if (norm_pl(n) > 1 - contrast/2) && (found_left)
            peak_right = n;
            break;
        end
    end
    
    FWHM = detune(peak_right) - detune(peak_left); % calculate FWHM
    
end
