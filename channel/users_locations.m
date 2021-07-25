
function [x_UE, y_UE] = users_locations(par)

min_angle_sep = par.MinAngleSep; % minumum angle seperation between users [deg]
broadside_BS_deg = 0; % broadside of BS antenna array  [deg]
% distance spread
if ~isfield(par, 'U0')
    par.U0 = par.U;
end

if par.U0 < par.U
    error('Total number of active users must be larger than the number of served users!');
end

d_min = par.dmin;
d_max = par.dmax;
d_UE = (d_max - d_min)*rand(par.U0,1) + d_min;

% ---------------------------------------------------------------------------
%--------------------- USERS DISTRIBUTION BEGINS HERE -----------------------
% GOAL: Generate par.U0 uniformly distributed random variables with minimum
% distance given by min_angle_sep

% Determine the angle of the first user

sector_edge1 = broadside_BS_deg-60;
sector_edge2 = broadside_BS_deg+60;

if par.MinAngleSep == 0
    aod_UE = unifrnd(sector_edge1,sector_edge2, par.U0,1);
else
    aod_UE = zeros(par.U0,1);
    p = [sector_edge1,sector_edge2];    % each row of p represents an unoccupied range of angles (now just one row!)
    aod_UE(1) = unifrnd(p(1), p(2));
    e1 = aod_UE(1) - min_angle_sep;
    e2 = aod_UE(1) + min_angle_sep;
    if e1 < sector_edge1
        p = [e2, p(2)];
    elseif e2 > sector_edge2
        p = [p(1), e1];
    else
        p = [p(1),e1; e2,p(2)];
    end
    
    % ----- Find the uniformly distributed random angle for each user (other than the first one)
    % such that the minimum angle seperation is met.
    for aod_i  = 2:par.U0
        % Find a uniformly distributed random angle for the next user
        % within the unoccupied ranges of angles
        unoccupied = zeros(size(p,1),1);
        for ii = 1:length(unoccupied)
            unoccupied(ii) = sum(p(1:ii,2) - p(1:ii,1));
        end
        raw_rand = unifrnd(0, unoccupied(end));
        aod_UE(aod_i) = p(1,1) + raw_rand;
        for ii = 1:length(unoccupied)-1
            if raw_rand > unoccupied(ii) && raw_rand < unoccupied(ii+1)
                aod_UE(aod_i) = p(ii+1,1) + raw_rand - unoccupied(ii);
            end
        end
        e1 = max(sector_edge1, aod_UE(aod_i) - min_angle_sep) + 0.0001;
        e2 = min(sector_edge2, aod_UE(aod_i) + min_angle_sep) - 0.0001;
        
        % -------- Determine each of e1 and e2 are within which range -----
        if e1 < p(1,1)
            e1p = -1;
        elseif e1 > p(1,1) && e1 < p(1,2)
            e1p = 1;
        end
        
        if e2 > p(end,2)
            e2p = -size(p,1);
        elseif e2 > p(end,1) && e2 < p(end,2)
            e2p = size(p,1);
        end
        for ii = 2:size(p,1)
            if e1 > p(ii,1) && e1 < p(ii,2)
                e1p = ii;
            elseif e1 > p(ii-1,2) && e1 < p(ii,1)
                e1p = -ii;
            end
        end
        
        
        for ii = 1:size(p,1)-1
            if e2 > p(ii,1) && e2 < p(ii,2)
                e2p = ii;
            elseif e2 > p(ii,2) && e2 < p(ii+1,1)
                e2p = -ii;
            end
        end
        % --------- e1 and e2 position determined! ---------------
        % --------- Now update the p matrix ---------------
        p_row = abs(e1p);
        if e1p > 0 && e2p > 0
            if e1p == 1
                p = [[p(1,1),e1];[e2,p(1,2)];p(2:end,:)];
            elseif e2p == size(p,1)
                p = [p(1:end-1,:);[p(end,1),e1];[e2,p(end,2)]];
            else
                p = [p(1:p_row-1,:);[p(p_row,1),e1];[e2,p(p_row,2)];p(p_row+1:end,:)];
            end
        elseif e1p < 0 && e2p < 0
            if e1p == -1
                p = p(2:end,:);
            elseif e1p == -size(p,1)
                p = p(1:end-1,:);
            else
                p = [p(1:p_row-1,:);p(p_row+1:end,:)];
            end
        elseif e1p < 0 && e2p > 0
            p(p_row,1) = e2;
        else
            p(p_row,2) = e1;
        end
    end
    
end
%---------------------- USERS DISTRIBUTION ENDS HERE ------------------------
% ---------------------------------------------------------------------------

if (length(aod_UE) > 1) && any(abs(diff(sort(aod_UE, 'ascend'))) < 0.99*min_angle_sep)
    error('Error in user locations: min angle separation not satisfied!');
end

coord_BS = [0, 0]; % BS coord.
coord_UE = ones(par.U0,1)*coord_BS + (d_UE*ones(1,2)).*[cosd(aod_UE), sind(aod_UE)]; % UE coord.

x_UE = coord_UE(:,1);
y_UE = coord_UE(:,2);

end
