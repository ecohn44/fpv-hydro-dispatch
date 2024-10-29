%% Lake Mead Data Import

% Total Capacity of Lake Mead (Ac-Ft)
mead_vol = [30167000; 28946000; 28667000; 26483000; 10230000; 4553000; 2547000];
meadlogvol = log(mead_vol);
mead_vol_m = mead_vol*1233.48; % convert Ac-Ft to cubic meters

% Elevation of Lake Mead (Ft)
mead_h = [1229.0; 1221.4; 1219.6; 1205.4; 1050.0; 950.0; 895.0];
meadlogh = log(mead_h);
mead_h_m = mead_h*0.3048; % convert ft to m

%% Lake Powell Data Import
df = readtable("data\Lake_Powell_2018_ElevAreaCap.csv");

% Capacity of Lake Powell (Ac-Ft)
pow_vol = df.Capacity_acrefeet;
powlogvol = log(pow_vol);
pow_vol_m = pow_vol*1233.48;

% Elevation of Lake Powell (Ft)
pow_h = df.Elevation_ft_NAVD88;
powlogh = log(pow_h);
pow_h_m = pow_h*0.3048;

%% Combined System 
pow_vol_subset_m = [pow_vol_m(883), pow_vol_m(1071), pow_vol_m(1379), pow_vol_m(1802), pow_vol_m(1821)]';
pow_h_subset_m = [pow_h_m(883), pow_h_m(1071), pow_h_m(1379), pow_h_m(1802), pow_h_m(1821)]';
mead_vol_subset_m = flip(mead_vol_m(3:7)); % last five values
mead_h_subset_m = flip(mead_h_m(3:7));

combo_vol_m = pow_vol_subset_m + mead_vol_subset_m;
combo_h_m = pow_h_subset_m + mead_h_subset_m;