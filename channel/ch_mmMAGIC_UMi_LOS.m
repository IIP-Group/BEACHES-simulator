%% Channel model setup and coefficient generation

function H = ch_mmMAGIC_UMi_LOS(par,x_UE,y_UE)

show = 0; % 0 = don't print status messages
% set up some basic parameters
if(show)
disp('set up scenario');
end
par.scenario = 'mmMAGIC_UMi_LOS'; % 'BERLIN_UMa_NLOS', 'Freespace', 'mmMAGIC_UMi_LOS', 'mmMAGIC_UMi_NLOS'
par.fc = 60e9; % carrier frequency [Hz]
par.BW = 40e6; % par.BW [Hz]
par.N = 128; % number of carriers
%par.B = 128; % number of BS antennas
%par.U = 8; % number of single-antenna UEs % number of user equipments

i_snapshot = 1; % used for frequency response: freq_response = c.fr(bandwidth,subcarriers,i_snapshot)

% Set up simulation parameters
s = qd_simulation_parameters;                           
s.show_progress_bars = false;                               
s.center_frequency = par.fc;                           
s.sample_density = 2.5; % 2.5 samples per half-wavelength
s.use_absolute_delays = false; % Include delay of the LOS path
s.use_spherical_waves = true; % Use spherical waves instead of plane waves?

lambda = s.wavelength;

%% Setup channel layout 

% Create new QuaDRiGa layout
l = qd_layout(s);    
if(show)
disp('generating positions');
end

% Define UE layout
l.no_rx = par.U;                                           
%%%% ???
for ii=1:l.no_rx
    % each user is static--> only 1 snapshot per track
    l.track(ii).no_snapshots = 1;
end
% place users at 1.5m
UE_pos_mat = [ x_UE' ; y_UE' ; repmat(1.5, 1, par.U) ];
% Write user locations into Quadriga layout object
l.rx_position = UE_pos_mat;
l.rx_array = qd_arrayant('omni'); % Omnidirectional UE antennas


% Set the scenario
l.set_scenario(par.scenario);


% Define BS antenna array
l.tx_array.generate('omni'); % so far, no need to use other antenna
l.tx_array.no_elements = par.B;
BS_pos_mat = zeros(3,l.tx_array.no_elements);
% positioning of array elements. units are in m. there is freedom to choose
% the position, multiply the spacing by lambda, etc.
for ii=1:l.tx_array.no_elements
    % disposing elements only on y-axis at half wavelength spacing
    BS_pos_mat(2,ii) = lambda/2*(ii-1) - lambda/2*(l.tx_array.no_elements-1)/2;    
end
l.tx_array.element_position = BS_pos_mat;
l.tx_position(3) = 25; % 25 m tx height

%% plot scenario
%disp('plotting scenario');
%l.visualize([],[],0);                                   % Plot the layout
%view(-33, 60);                                          % Enable 3D view

%% Generate channel coefficients

% Instead of creating a channel builder, can just call l.get_channels from 
% the layout object, which will do the same thing as below since we do not
% change any of the scenario (Berlin) parameters.
%{
p = l.init_builder; % Initialize builder
p.gen_ssf_parameters;                                  
p.gen_lsf_parameters;
%}

if(show)
disp('generate channels');
end
c = l.get_channels; % Generate channels

%% Calculating the frequency response for this channel
if(show)
disp('genereate frequency response');
end

% Writing the frequency response in the desired way: a matrix corresponding
% to the user grid, each point with a vector with the number of antennas
HN = zeros(par.B,par.U,par.N);
for uu=1:par.U
    HN(:,uu,:) = c(uu).fr(par.BW,par.N,i_snapshot); % Using built-in method to evaluate the frequency response from the channel coefficients
end

if(show)
disp('extract channel')
end

% just extract the first subcarrier
H = HN(:,:,1);

end
