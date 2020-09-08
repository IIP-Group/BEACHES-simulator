% -----------------------------------------------------
% -- LoS MIMO channel with some minumum angle seperation between users 
% -- 2018 (c) studer@cornell.edu, sm2675@cornell.edu
% -- Last update: 12/03/2018
% -----------------------------------------------------

function [H_swm, H_pwm] = ch_LoS(par, x_UE, y_UE)

    c = 3e8; % speed of light [m/s]
    lambda =  c / 2.4e9; % carrier wavelength [m]
    delta = .5; % antenna spacing
    broadside_BS_deg = 0; % broadside of BS antenna array  [deg]

    d_avg = mean(sqrt(x_UE.^2 + y_UE.^2));

    coord_BS = [0, 0]; % BS coord.
    coord_UE(:,1) = x_UE;
    coord_UE(:,2) = y_UE;
    
    d_ant_BS = delta * lambda; % distance between BS antenna elements [m]
    d_array_BS = (par.B - 1) * d_ant_BS; % size of BS antenna array [m]

    % array rotation
    Omega_BS_deg = wrapTo360(90 - broadside_BS_deg); % BS array rotation [deg]
    Omega_BS_rad = pi/180 * Omega_BS_deg; % BS array rotation [rad]
    Omega_UE_deg = wrapTo360(unifrnd(-180, 180, par.U, 1) - 90); % UE array rotation [deg]
    Omega_UE_rad = pi/180 * Omega_UE_deg; % UE array rotation [rad]

    % antenna elem. coordinates
    x_BS = coord_BS(1) + d_ant_BS*((1:par.B) - (par.B+1)/2)*cos(pi-Omega_BS_rad);
    y_BS = coord_BS(2) + d_ant_BS*((1:par.B) - (par.B+1)/2)*sin(pi-Omega_BS_rad);
    

    % LoS AoD
    theta_LoS_BS_rad = Omega_BS_rad - pi/2 + atan2((coord_UE(:,2) - coord_BS(2)),(coord_UE(:,1) - coord_BS(1)));
    theta_LoS_BS_deg = 180/pi * theta_LoS_BS_rad; 
    theta_LoS_UE = pi - Omega_UE_rad - atan2((coord_UE(:,1) - coord_BS(1)),(coord_UE(:,2) - coord_BS(2)));
    theta_LoS_UE_deg = 180/pi * theta_LoS_UE; 

    % coordinates
    xx_BS = ones(par.U,1)*x_BS; yy_BS = ones(par.U,1)*y_BS;
    xx_UE = x_UE*ones(1,par.B); yy_UE = y_UE*ones(1,par.B);
    
    % reference distance
    d_ref = sqrt((xx_BS - xx_UE).^2 + (yy_BS - yy_UE).^2);

    % angles
    theta_BS = Omega_BS_rad - pi/2 + atan2((yy_UE - yy_BS),(xx_UE - xx_BS));
    theta_UE = pi - Omega_UE_rad*ones(1,par.B) - atan2((xx_UE - xx_BS),(yy_UE - yy_BS));
    
    % distances between UE and BS antenna elements
    dd_ant_BS = d_ant_BS*ones(par.U,1)*((1:par.B)-1); tt_BS = theta_BS(:,1)*ones(1,par.B);
    
    % distance according to PWM model
    d_pwm = d_ref(:,1)*ones(1,par.B) - dd_ant_BS.*sin(tt_BS);
    
    % distance according to SWM model
    d_swm = sqrt((d_ref(:,1)*ones(1,par.B) - dd_ant_BS.*sin(tt_BS)).^2 + (dd_ant_BS.*cos(tt_BS)).^2);

    if sum(sum(abs(d_ref - d_swm) < 0.001./d_ref)) ~= par.B*par.U
       warning('something is wrong here...'); 
    end

    % channel matrix
    H_pwm = d_avg./d_pwm .* exp(-1i*2*pi*d_pwm/lambda);
    H_swm = d_avg./d_swm .* exp(-1i*2*pi*d_swm/lambda);
    H_pwm = H_pwm.';
    H_swm = H_swm.';
%     H_pwm  = H_pwm'/sqrt(mean(abs(H_pwm(:)).^2));
%     H_swm  = H_swm'/sqrt(mean(abs(H_swm(:)).^2));
end

function lon = wrapTo360(lon)

    positiveInput = (lon > 0);
    lon = mod(lon, 360);
    lon((lon == 0) & positiveInput) = 360;

end