function PL = single_nv_odmr(gamma, w1, rabi_freq, optical_power, d)
    %% INPUT PARAMETERS
    % gamma         = 1/T1 (Hz)
    % w1            = frequency splitting between 
    % rabi_freq     = measured Rabi oscillation frequency (rads/sec)
    % optical_power = excitation rate of 532nm laser
    % d             = microwave detuning (Hz)
    
    %% Modeling 5 level N-V ODMR with hyperfine and power broadening effects % % % 
    %
    %
    %                        4 ----  
    %     3 ----     
    %
    %              5 ----
    %
    %                        2 ---- 
    %     1 ---- 
    %   
    % 1: ms = 0   ground     state  |3A2,0>
    % 2: ms = -1  ground     state  |3A2,-1>
    % 3: ms = 0   excited    state  |3E,0>
    % 4: ms = -1  excited    state  |3E,-1>
    % 5: shelf    metastable state  |1E> 

    %% EMPIRICAL CONSTANTS
    
    % w1     = 2.87 * 10^9;         % energy splitting (Hz)
    hbar   = 1; % 1.0545 * 10^(-34);     % (m^2 kg/s)
    WR     = rabi_freq;             % spin Rabi frequency (placeholder value) (Hz)
    G10    = 1 / (1 * 10^(-3));     % population decay from |1> to |0> %ranges from 3 milliseconds to 4 ms according to Eisuke's measurements on the Itoh sample (Hz)
    G01    = G10;                   % population fluctuation from |0> to |1> (G01 = G10 e^(-dE/kbT), where dE/kbT ? 0) (Hz)
    G20    = 1 / (10 * 10^(-9));    % population decay from |2> to |0> (Hz)
    P02    = optical_power;         % |0> to |2> pump rate, proportional to optical power (arbitrary)
    P13    = P02;                   % |1> to |3> pump rate, proportional to optical power (arbitrary)
    G31    = G20;                   % population decay from |3> to |1> (Hz)
    G34    = G20;                   % population decay from |3> to |4> (Hz) (on the order of G20)
    G40    = 1 / (300 * 10^(-9));   % population decay from |4> to |0> (Hz)
    % gamma  = 1 / (1 * 10^(-6));   % dephasing rate (1/T2 process) (Hz)

    
    %% SOLVE FOR STEADY STATE SOLUTION TO RATE EQUATIONS
    
    % formed by (i/hbar)[H0,rho] + L(rho) = 0 (steady state)
        % |00>                                                                 |11>       |22>     |33>   total pop.
    M = [-P02-G01           -(1i/hbar)*WR/2         (1i/hbar)*WR/2             G10         G20      0     G40;  % |0>
        -(1i/hbar)*WR/2    -gamma+(1i/hbar)*(d-w1)          0             (1i/hbar)*WR/2   0        0      0;
        (1i/hbar)*WR/2              0            -gamma-(1i/hbar)*(d-w1) -(1i/hbar)*WR/2   0        0      0; 
        G01                  (1i/hbar)*WR/2        -(1i/hbar)*WR/2          -P13-G10       0       G31     0;   % |1>
        P02                         0                       0                   0        -G20       0      0;   % |2>
        0                           0                       0                  P13         0     -G34-G31  0;   % |3>
        1                           0                       0                   1          1        1      1];  % total pop.

    % steady state vector
    steadystate = [0,0,0,0,0,0,1]';

    % solve for rho when drhodt == 0
    rho = M \ steadystate;

    % PL intensity (observable)
    PL = real(G20*rho(5) + G31*rho(6));
    
end