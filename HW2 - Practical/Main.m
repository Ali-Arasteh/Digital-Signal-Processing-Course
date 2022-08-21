%% ==== Loading Data ====

disp('Loading Signals...')
[hdr, records] = edfread('ST7011J0-PSG.edf');
disp('Done');

%% ==== Filtering ====
 
Fpz = records(1,1:end);
Oz  = records(2,1:end);

% your turn :)