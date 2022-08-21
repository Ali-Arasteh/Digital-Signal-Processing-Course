%% ==== Loading Data ====

disp('Loading Signals...')
[hdr, records] = edfread('ST7011J0-PSG.edf');
disp('Done');

%% ==== Filtering ====

disp('Filtering Frequency Bands...')
h_alpha  = BPF (1001, 8  , 15, 100);
h_beta   = BPF (1001, 16 , 31, 100);
h_theta  = BPF (1001, 4  ,  7, 100);
h_delta  = BPF (1001, 0.5,  4, 100);
 
Fpz = records(1,1:end);
Oz  = records(2,1:end);

Fpz_alpha = filter(h_alpha, 1, Fpz);
Fpz_beta  = filter(h_beta , 1, Fpz);
Fpz_theta = filter(h_theta, 1, Fpz);
Fpz_delta = filter(h_delta, 1, Fpz);

Oz_alpha  = filter(h_alpha, 1, Oz);
Oz_beta   = filter(h_beta , 1, Oz);
Oz_theta  = filter(h_theta, 1, Oz);
Oz_delta  = filter(h_delta, 1, Oz);

disp('Done')

%% ==== Ploting ====
% your turn!