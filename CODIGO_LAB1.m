%%%%LABORATORY SESSION 1%%%%
clc;clear;close all;

% Orbit parameters----------------------------------

% Orbit A
a_a = 26556e3; %semi-mayor axis
e_a = 0.70399; %eccentricity
i_a = 63.4; %inclination
raan_a = 90; %right ascension of the ascending node
aop_a = 270; %argument of periapsis
ta_a = 0; %true anomaly

% Orbit B
a_b = 7039.5e3; %semi-mayor axis
e_b = 0.0001; %eccentricity
i_b = 98.074; %inclination
raan_b = 287; %right ascension of the ascending node
aop_b = 260.31; %argument of periapsis
ta_b = 0; %true anomaly

% Orbit C
a_c = 7202.5e3; %semi-mayor axis
e_c = 0.02; %egocentricity
i_c = 42.2; %inclination
raan_c = 121; %right ascension of the ascending node
aop_c = 15; %argument of periapsis
ta_c = 0; %true anomaly

% Orbit D
a_d = 6633.137e3; %semi-mayor axis
e_d = 0.0009826; %eccentricity
i_d = 96.5717; %inclination
raan_d = 344.5256; %right ascension of the ascending node
aop_d = 296.2811; %argument of periapsis
ta_d = 0; %true anomaly

names = ["SAT_A","SAT_B","SAT_C","SAT_D"];
%       a [m]        e          i [deg]   RAAN [deg]  argp [deg]  nu [deg]
elems = [26556e3     0.70399    63.4      90          270         0; % A:1
        7039.5e3    0.0001     98.074    287         260.31      0;  % B:2
        7202.5e3    0.02       42.2      121         15          0;  % C:3
        6633.137e3  0.0009826  96.5717   344.5256    296.2811    0]; % D:4

% T0
start_time = datetime(2013,11,10,23,03,03.95,'TimeZone','UTC');

% Global parameters
mu = 3.986004418e14;        
Re = 6371137;                 % [km]
J2 = 1.08262668e-3;
AU  = 149597870.7;            % astronomical unit [km]

%Satellite properties 
mass_sat  = 900.0;      % Satellite mass [kg]
area_drag = 1.1;        % Cross-sectional area S [m^2]
Cd        = 3.0;       % Drag coefficient [-]

%% Task 1(a)
clc;
close all;

% Select end date and sample time 
stop_time = start_time+days(10); %propagate for 10 day
sample_time=120; %Sample every 60 seconds

% Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

% Numerical propagator:
numericalPropagator(sc,...
    'IncludeAtmosDrag', false, ...
    'GravitationalPotentialModel', "oblate-ellipsoid", ... %
    'SphericalHarmonicModel','EGM2008',...
    'ODESolver', "ode45", ... 
    'ODESet', odeset('RelTol', 1e-8, 'AbsTol', 1e-8));

% Satellites 
% satA = satellite(sc,a_a,e_a,i_a,raan_a,aop_a,ta_a,"Name","SAT_A", ...
%     "OrbitPropagator","two-body-keplerian");
satB = satellite(sc,a_b,e_b,i_b,raan_b,aop_b,ta_b,"Name","SAT_B", ...
    "OrbitPropagator","two-body-keplerian");
% satC = satellite(sc,a_c,e_c,i_c,raan_c,aop_c,ta_c,"Name","SAT_C", ...
%     "OrbitPropagator","numerical");
% satD = satellite(sc,a_d,e_d,i_d,raan_d,aop_d,ta_d,"Name","SAT_D", ...
%     "OrbitPropagator","numerical");

% Visualize orbits
%satelliteScenarioViewer(sc,"Dimension","3D");
% satelliteScenarioViewer(sc,"Dimension","2D");

% Ground Track
% groundTrack(satA,LeadTime=86400*1);
% groundTrack(satB,LeadTime=86400*1);
% groundTrack(satC,LeadTime=86400*10);
% groundTrack(satD,LeadTime=86400*1);

% Check if sunsynchronous orbits
% isSunSynchronous(elems(2,:), names(2), 0.01, start_time, Re, mu, J2);
% isSunSynchronous(a_b, e_b, i_b, "Orbit B", 0.01);
% isSunSynchronous(a_c, e_c, i_c, "Orbit C", 0.01);
% isSunSynchronous(a_d, e_d, i_d, "Orbit D", 0.01);

for i=1:4
    kep = elems(i,:);

    % Orbital characteristics
    [ha, hp, T] = getOrbitCharacteristics(kep, Re, mu);
    fprintf("\n===== Orbit %s =====\n", names(i));
    fprintf('Apogee Altitude : %.1f km\n', ha);
    fprintf('Perige Altitude : %.1f km\n', hp);
    fprintf('Apogee Altitude : %.1f hours\n', T/3600);

    % Sun-synchronous?
    epoch = datetime(2013,11,10,23,03,03.95);
    isSunSynchronous(kep, names(i), 0.01, epoch, Re, mu, J2);
end



%% Task 1(b)
clc;
close all;

% Time propagation settings
total_duration = days(365);
stop_time = start_time + total_duration;
 
% Sample time per orbit [A B C D] in s 
sample_times = [30 30 30 30];
% Maximum samples per chunk. 
max_samples = 20000;         
 
% Select orbit/orbits to run. Recommended one orbit at a time
orbitsToRun = 1; 
 
% Ground station data
gs_lat = 40;             
gs_lon = -3.7033;        
min_elev = 20;             
 
% Preallocate results 
nOrb = numel(orbitsToRun);
nContacts = zeros(nOrb,1);
avgWin = nan(nOrb,1);   % [min]
minWin = nan(nOrb,1);   % [min]
maxWin = nan(nOrb,1);   % [min]
timePerDay = nan(nOrb,1);   % [min/day] 
maxGap = nan(nOrb,1);   % [h] longest time without contact
 
% Main loop: orbits -> chunks
for k = 1:nOrb
    idx = orbitsToRun(k);
    sample_time  = sample_times(idx);
    chunk_length = seconds(max_samples*sample_time);
    fprintf("\n===== Orbit %s (step %g s, chunks of %.2f days) =====\n", ...
        names(idx), sample_time, days(chunk_length));
 
    kep = elems(idx,:);    % [a e i RAAN argp nu] at the start of the current chunk
    t0  = start_time;      % start of the current chunk
 
    % Access intervals of the whole year, stored as seconds since epoch
    allStart = zeros(0,1);
    allEnd   = zeros(0,1);
 
    chunk = 0;
    while t0 < stop_time
        chunk = chunk + 1;
        t1 = min(t0 + chunk_length, stop_time);   % last chunk may be shorter
 
        % Satellite scenario for this chunk
        sc = satelliteScenario(t0, t1, sample_time);
        numericalPropagator(sc, ...
            'IncludeAtmosDrag', false, ...
            'GravitationalPotentialModel', "oblate-ellipsoid", ... % J2
            'ODESolver', "ode45", ...
            'ODESet', odeset('RelTol',1e-8,'AbsTol',1e-8));
 
        % Satellite initialized with the state at t0
        sat = satellite(sc, kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), ...
            "Name", names(idx), "OrbitPropagator", "numerical");
 
        % Ground station and access
        gs = groundStation(sc, gs_lat, gs_lon, ...
            "MaskElevationAngle", min_elev, "Name", "GS");
        ac     = access(sat, gs);
        intvls = accessIntervals(ac);
 
        % Check: elevation at the start of thefirst pass should be aroud 20 deg
        if chunk == 1 && height(intvls) > 0
            [~, el] = aer(gs, sat, intvls.StartTime(1));
            fprintf("Check: elevation at first pass start = %.2f deg\n", el);
        end
 
        %  Save the relevant info and merge pass if cut in two by chunk boundary
        s = intvls.StartTime;
        s.TimeZone = 'UTC';
        e = intvls.EndTime;
        e.TimeZone = 'UTC';
        newStart = seconds(s - start_time);
        newEnd   = seconds(e - start_time);
        chunkStartSec = seconds(t0 - start_time);
        [allStart, allEnd] = appendIntervals(allStart, allEnd, ...
            newStart, newEnd, chunkStartSec, sample_time/2);
 
        % Prepare next chunk propagation
        [r_end, v_end] = states(sat, t1, "CoordinateFrame", "inertial");
        [a, ecc, incl, RAAN, argp, nu] = ijk2keplerian(r_end, v_end);
        kep = [a ecc incl RAAN argp nu];
 
        fprintf("Chunk %2d: %s -> %s ; %3d passes\n", chunk, ...
            string(t0,'dd-MMM-yyyy'), string(t1,'dd-MMM-yyyy'), height(intvls));
 
        t0 = t1;
        clear sc sat gs ac intvls   % free memory before the next chunk
    end
 
    % Statistics for this orbit
    dur_min = (allEnd - allStart)/60;                 % window durations [min]
    nContacts(k)  = numel(dur_min);
    avgWin(k)     = mean(dur_min);
    minWin(k)     = min(dur_min);
    maxWin(k)     = max(dur_min);
    timePerDay(k) = sum(dur_min)/days(total_duration);
    if numel(allStart) > 1
        maxGap(k) = max(allStart(2:end) - allEnd(1:end-1))/3600;
    end
 
end
 
% Results
contactsPerDay = nContacts/days(total_duration);
Results = table(names(orbitsToRun).', nContacts, contactsPerDay, avgWin, ...
    minWin, maxWin, timePerDay, maxGap, 'VariableNames', ...
    {'Orbit','Contacts_year','Contacts_per_day','AvgWindow_min', ...
     'MinWindow_min','MaxWindow_min','ContactTime_min_per_day','MaxGap_h'});
disp(Results)

%% Task 1(c) 
clc; close all;

% Settings
stop_time  = start_time + days(30);      

SelectOrbit = 4;   % 1 = A, 2 = B, 3 = C, 4 = D
sample_times = [30 15 15 10];         % step in s for orbits A, B, C, D
sample_time  = sample_times(SelectOrbit);
% Maximum samples per chunk. 
max_samples = 20000;
% Chunk length
chunk_length = seconds(max_samples*sample_time);  

% State at the start of each chunk
kep = elems(SelectOrbit,:);               

eq_lon = -180:22.5:135;          % [deg] longitudes of the equator points
nPts   = numel(eq_lon);

% Access intervals of each point (seconds since start_time)
intS = cell(1, nPts);
intE = cell(1, nPts);

% Propagation in chunks
t0 = start_time;
while t0 < stop_time
    t1 = min(t0 + chunk_length, stop_time);

    % Scenario and propagator for this chunk
    sc = satelliteScenario(t0, t1, sample_time);
    numericalPropagator(sc, 'IncludeAtmosDrag', false, ...
        'GravitationalPotentialModel', "oblate-ellipsoid", ...
        'ODESolver', "ode45", 'ODESet', odeset('RelTol',1e-8,'AbsTol',1e-8));
    sat = satellite(sc, kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), ...
        "OrbitPropagator", "numerical");

    % Nadir-pointing camera
    g = gimbal(sat);
    pointAt(g, "nadir");
    cam = conicalSensor(g, "MaxViewAngle", 40);   

    % Fake ground stations on the equator and camera access to each one
    ac = cell(1, nPts);
    for p = 1:nPts
        pt    = groundStation(sc, 0, eq_lon(p));
        ac{p} = access(cam, pt);
    end

    % Store the intervals of this chunk
    for p = 1:nPts
        iv = accessIntervals(ac{p});
        S = iv.StartTime;  S.TimeZone = 'UTC';
        E = iv.EndTime;    E.TimeZone = 'UTC';
        [intS{p}, intE{p}] = appendIntervals(intS{p}, intE{p}, ...
            seconds(S - start_time), seconds(E - start_time), ...
            seconds(t0 - start_time), sample_time/2);
    end

    % Prepare next chunk propagation
    [r_end, v_end] = states(sat, t1, "CoordinateFrame", "inertial");
    [a, ecc, incl, RAAN, argp, nu] = ijk2keplerian(r_end, v_end);
    kep = [a ecc incl RAAN argp nu];

    fprintf("Done until %s\n", string(t1, 'dd-MMM-yyyy HH:mm'));
    t0 = t1;
    clear sc sat g cam ac pt iv
end

% Revisit time
maxRevisit_h = zeros(nPts, 1);
nAccess      = zeros(nPts, 1);
for p = 1:nPts
    nAccess(p) = numel(intS{p});
    if nAccess(p) >= 2
        gaps = intS{p}(2:end) - intE{p}(1:end-1);   % [s] end of one access, start of next
        maxRevisit_h(p) = max(gaps)/3600;
    else
        maxRevisit_h(p) = Inf;                       % seen fewer than 2 times
    end
end

disp(table(eq_lon.', nAccess, maxRevisit_h, ...
    'VariableNames', {'Lon_deg','nAccess','MaxRevisit_h'}))
fprintf("Maximum revisit at the equator = %.2f h (%.2f days)\n", ...
    max(maxRevisit_h), max(maxRevisit_h)/24);

%% Task 1(d)
clc;close all;

%Select end date and sample time 
stop_time = start_time+days(10); %propagate for 10 day
sample_time=60; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

% Numerical propagator:

numericalPropagator(sc,...
    'IncludeAtmosDrag', false, ...
    'GravitationalPotentialModel', "oblate-ellipsoid", ... %
    'SphericalHarmonicModel','EGM2008',...
    'ODESolver', "ode45", ... 
    'ODESet', odeset('RelTol', 1e-8, 'AbsTol', 1e-8)); 

% Orbit selection=0 -> Orbit A. Any thing else Orbit C
OrbitSelection = 0; 
if OrbitSelection == 0
    sat = satellite(sc,a_a,e_a,i_a,raan_a,aop_a,ta_a,"Name","SAT_A", ...
        "OrbitPropagator","two-body-keplerian");
    sat_J2 = satellite(sc,a_a,e_a,i_a,raan_a,aop_a,ta_a,"Name","SAT_A_J2", ...
        "OrbitPropagator","numerical");
    PlotOrbit = "A";
else
    sat = satellite(sc,a_c,e_c,i_c,raan_c,aop_c,ta_c,"Name","SAT_C", ...
        "OrbitPropagator","two-body-keplerian");
    sat_J2 = satellite(sc,a_c,e_c,i_c,raan_c,aop_c,ta_c,"Name","SAT_C_J2", ...
        "OrbitPropagator","numerical");
    PlotOrbit = "C";
end


% Plot orbital elements wit and without J2---------------------------------
% Get orbital arc
[r_noJ2, v_noJ2, t] = states(sat);      % 3xN [m], 3xN [m/s], 1xN datetime
[r_J2,   v_J2]      = states(sat_J2);

% Convert to keplerian elements
[a_noJ2,ecc_noJ2,incl_noJ2, raan_noJ2,argp_noJ2, nu_noJ2,~,~] = ...
    ijk2keplerian(r_noJ2, v_noJ2, GravitationalParameter=mu);
[a_J2,ecc_J2,incl_J2, raan_J2, argp_J2, nu_J2,~,~] = ...
    ijk2keplerian(r_J2, v_J2, GravitationalParameter=mu);

% Plot orbital elements vs time
tdays = days(t - t(1));                       % time since epoch [days]

figure('Name', "Orbit " + PlotOrbit + " - orbital elements");
tl = tiledlayout(3,2,'TileSpacing','compact');

nexttile
plot(tdays, a_noJ2/1e3, tdays, a_J2/1e3)
ylabel('a [km]'); grid on
legend('Without J_2','With J_2','Location','best')

nexttile
plot(tdays, ecc_noJ2, tdays, ecc_J2)
ylabel('Eccentricity'); grid on
legend('Without J_2','With J_2','Location','best')

nexttile
plot(tdays, incl_noJ2, tdays, incl_J2)
ylabel('Inclination [deg]'); grid on
legend('Without J_2','With J_2','Location','best')

nexttile
plot(tdays, raan_noJ2, tdays, raan_J2)
ylabel('RAAN [deg]'); grid on
legend('Without J_2','With J_2','Location','best')

nexttile
plot(tdays, argp_noJ2, tdays, argp_J2)
ylabel('Argument of perigee [deg]'); grid on
legend('Without J_2','With J_2','Location','best')

nexttile
plot(tdays, nu_noJ2, tdays, nu_J2)
ylabel('True Anomaly [deg]'); grid on
legend('Without J_2','With J_2','Location','best')

xlabel(tl, 'Time since epoch [days]')
title(tl, "Orbit " + PlotOrbit + ": osculating orbital elements")

% Plot ground track--------------------------------------------------------
% Get lat and long
lla_noJ2 = states(sat, "CoordinateFrame","geographic");   
lat_noJ2 = lla_noJ2(1,:);                                 
lon_noJ2 = lla_noJ2(2,:);                                

lla_J2 = states(sat_J2, "CoordinateFrame","geographic");  
lat_J2 = lla_J2(1,:);                                 
lon_J2 = lla_J2(2,:);                                

% Insert NaN wherever longitude jumps by more than 180 deg between samples
[latP_noJ2, lonP_noJ2] = breakAtDateline(lat_noJ2, lon_noJ2);
[latP_J2, lonP_J2] = breakAtDateline(lat_J2, lon_J2);

figure('Name', "Orbit " + PlotOrbit + " - ground track");

% Coastlines: 0 m contour of the base-MATLAB topography data
load topo                                         
topo = [topo(:,181:360), topo(:,1:180)];          
contour(-180:179, -89:90, topo, [0 0], ...
    'Color',[0.6 0.6 0.6], 'HandleVisibility','off');
hold on

% Ground tracks plot
plot(lonP_noJ2, latP_noJ2, '-',  'LineWidth',1.0, 'DisplayName','Without J_2')
plot(lonP_J2,   latP_J2,   '--', 'LineWidth',1.0, 'DisplayName','With J_2')
plot(lon_noJ2(1), lat_noJ2(1), 'ko', 'MarkerFaceColor','k', 'DisplayName','Start')
plot(-3.7033, 40, 'k^', 'MarkerFaceColor','y', 'DisplayName','Ground station')

axis([-180 180 -90 90]); daspect([1 1 1])
xticks(-180:60:180); yticks(-90:30:90); grid on
xlabel('Longitude [deg]'); ylabel('Latitude [deg]')
title("Orbit " + PlotOrbit + ": ground track with and without J_2")
legend('Location','southoutside','Orientation','horizontal')
hold off

%% Task 1(e)
clc;close all;

%Select end date and sample time 
stop_time = start_time+days(10); 
sample_time=300; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

% Numerical propagator and physical properties:

numericalPropagator(sc,...
    'IncludeAtmosDrag', true, ...
    'GravitationalPotentialModel', "point-mass", ... %.oblate-ellipsoid
    'ODESolver', "ode45", ... 
    'ODESet', odeset('RelTol', 1e-7, 'AbsTol', 1e-7)); 
satD = satellite(sc,a_d,e_d,i_d,raan_d,aop_d,ta_d,"Name","SAT_D", ...
    "OrbitPropagator","numerical");
% satB = satellite(sc,a_b,e_b,i_b,raan_b,aop_b,ta_b,"Name","SAT_B", ...
%     "OrbitPropagator","numerical");
physProps = physicalProperties (satD,'Mass', mass_sat, ...
    'DragArea', area_drag, ...
    'DragCoefficient', Cd);


[r, v, t] = states(satD);
tsec  = seconds(t - t(1));              % time since epoch [s]
tsec  = tsec(:);
[a, ~, ~, ~, ~, ~, ~, ~] = ijk2keplerian(r, v, GravitationalParameter=mu);
a = a(:);                               

%-Secular decay rate of a  as a linear fit of the first two days
fit_window = 5*86400;                   
idx  = tsec <= fit_window;
p    = polyfit(tsec(idx), a(idx), 1);
adot = p(1);                            % [m/s]

fprintf('Decay rate of a : %.3f km/day\n', adot*86400/1e3);

% DeltaV needed to keep da/dt = 0
n0        = sqrt(mu/a_d^3);             
T_orbit   = 2*pi/n0;                    
dVdt      = 0.5*n0*abs(adot);           % [m/s^2]

T_mission = 8*365.25*86400;             % 8 years [s]
dV_orbit  = dVdt*T_orbit;
dV_day    = dVdt*86400;
dV_year   = dVdt*365.25*86400;
dV_total  = dVdt*T_mission;

fprintf('DeltaV per orbit : %.4f m/s\n', dV_orbit);
fprintf('DeltaV per day   : %.3f m/s\n', dV_day);
fprintf('DeltaV per year  : %.1f m/s\n', dV_year);
fprintf('DeltaV 8 years   : %.1f m/s\n', dV_total);

% Analytical cross-check with the density table                   
h0    = a_d - Re;
rho   = 6.16e-11;  % Density at 240km with mean solar activity
v0    = sqrt(mu/a_d);
B     = Cd*area_drag/mass_sat;          % ballistic coefficient
acc_drag      = 0.5*rho*v0^2*B;           % [m/s^2]
dV_total_an = acc_drag*T_mission;
da_rev      = -2*pi*rho*a_d^2*B;        % decay of a per revolution [m]

fprintf('\nAnalytical check');
fprintf('Decay per orbit  : %.1f m\n',   da_rev);
fprintf('DeltaV 8 years   : %.1f m/s\n', dV_total_an);

% Plot
figure('Name', "Orbit D - decay and fit");
plot(tsec/86400, a/1e3, 'b'); hold on
plot(tsec(idx)/86400, polyval(p, tsec(idx))/1e3, 'r--', 'LineWidth', 1.5)
xlabel('Time since epoch [days]'); ylabel('a [km]'); grid on
legend('Propagated a', 'Linear fit (secular rate)')
title('Orbit D: semi-major axis decay with drag')

%% Task 3(a); Chosen Orbit => B
% Report values obtained with 30 second sample time over a year. Task takes 
% approx 40 mins to run at that sample time. It is recommended to use a 
% higher sample time or lower stop time just to check if the program works
% if task3a_orbitB.mat table with the propagated values won't be used
clc; close all;

% Settings
idx = 2;                      % orbit B
sample_time  = 30;            % seconds
max_samples  = 20000;
chunk_length = seconds(max_samples*sample_time);
stop_time    = start_time + days(365);

kep   = elems(idx,:);
t0    = start_time;
chunk = 0;

if isfile('task3a_orbitB.mat')
    load('task3a_orbitB.mat');
else
% Chunk Propagation 
T_all    = [];          % Time of every sample in seconds
lit_all  = false(1,0);  % Eclipse logical of ev. sample, True = eclipse
dAU_all  = [];          % earth-Sun distance ev. sample [AU]
sinB_all = [];          % beta value of every sample
g_all = [];             % Shadow function
while t0 < stop_time
    chunk = chunk + 1;
    t1 = min(t0 + chunk_length, stop_time);

    % Propagation 
    sc = satelliteScenario(t0, t1, sample_time);
    numericalPropagator(sc, 'IncludeAtmosDrag', false, ...
        'GravitationalPotentialModel', "oblate-ellipsoid", ...
        'ODESolver', "ode45", 'ODESet', odeset('RelTol',1e-8,'AbsTol',1e-8));
    sat = satellite(sc, kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), ...
        "Name", names(idx), "OrbitPropagator", "numerical");

    % Satellite positions AND the sample times they belong to
    [r_sat, v_sat, t] = states(sat, "CoordinateFrame", "inertial"); % [m]
    t.TimeZone = 'UTC';
    r_km = r_sat / 1000; 
    tS = seconds(t - start_time);           % seconds since epoch

    %  Sun position at the same instants
    jd = juliandate(t);                   
    % Return position of Sun relative to earth centre (points to sun)
    r_sun = planetEphemeris(jd(:), 'Earth', 'Sun').';   % [km] 
    d_sun = vecnorm(r_sun);                             % distance from sun
    s_hat = r_sun ./ d_sun;                             % 3xN unit vectors
    d_AU = d_sun / AU;                                 % 1xN distance in AU

    % Cylindrical shadow eclipse model
    p = sum(r_km .* s_hat, 1);            % r·s_hat for every column[km]
    d_perp = vecnorm(r_km - p .* s_hat);  % 1xN distance to Earth-Sun line 
                                          %(perpendicular distance to cylinder axis) 
    % If satellite behind terminator line (p<0) and distance to earth sun line 
    % < Earth radius -> Eclipse; 1XN logical variable. true = eclipse
    inEcl  = (p < 0) & (d_perp < Re/1000);   

    % Refine entry/exit times by linear interpolation of g
    % shadow function: >0 outside cylinder, <0 inside
    g = d_perp - Re/1000; 

    % orbit normal and beta angle at every sample
    h_hat = cross(r_sat, v_sat, 1);          % column r x v
    h_hat = h_hat ./ vecnorm(h_hat);         % unit orbit normal
    sinB  = sum(h_hat .* s_hat, 1);          % sin(beta) = h_hat · s_hat

    % Skip first sample of every chunk after the first
    if chunk == 1
        keep = 1:numel(tS); 
    else
        keep = 2:numel(tS);
    end

    T_all    = [T_all,    tS(keep)];        
    lit_all  = [lit_all,  ~inEcl(keep)];    
    dAU_all  = [dAU_all,  d_AU(keep)];      
    sinB_all = [sinB_all, sinB(keep)];      
    g_all    = [g_all,    g(keep)];

    % Prepare next chunk propagation 
    [r_end, v_end] = states(sat, t1, "CoordinateFrame", "inertial");
    [a, ecc, incl, RAAN, argp, nu] = ijk2keplerian(r_end, v_end);
    kep = [a ecc incl RAAN argp nu];

    fprintf("Done until %s\n", string(t1, 'dd-MMM-yyyy HH:mm'));
    t0  = t1;
    clear sc sat
end
save('task3a_orbitB.mat', 'T_all', 'lit_all', 'dAU_all', 'sinB_all', 'g_all', 'sample_time');
end

% Eclipse interval interpolation over the full mission time (year)
dt  = sample_time;
N   = numel(T_all);
edges  = diff([false, ~lit_all, false]);        % +1 entry, -1 exit
iStart = find(edges ==  1);
iEnd   = find(edges == -1) - 1;

eclStart = T_all(iStart);                        % defaults: sample times
eclEnd   = T_all(iEnd);

m = iStart > 1;  k = iStart(m);                  % refine entries
eclStart(m) = T_all(k-1) + dt * g_all(k-1) ./ (g_all(k-1) - g_all(k));
m = iEnd < N;    k = iEnd(m);                    % refine exits
eclEnd(m)   = T_all(k)   + dt * g_all(k)   ./ (g_all(k)   - g_all(k+1));

full = (iStart > 1) & (iEnd < N);                % drop eclipses cut by start/end of run
eclStart = eclStart(full);
eclEnd   = eclEnd(full);
eclDur   = (eclEnd - eclStart) / 60;             % [min], no +dt correction needed

fprintf("Eclipses: %d | mean %.1f | min %.1f | max %.1f min\n", ...
    numel(eclDur), mean(eclDur), min(eclDur), max(eclDur));
figure; plot(eclStart/86400, eclDur, '.');
xlabel('Days since epoch'); ylabel('Eclipse duration [min]'); grid on;

% Power computation
P0   = 2000;                          % [W] at 1 AU if normal to sun
flux = lit_all ./ dAU_all.^2;         % 0 in eclipse, 1/d^2 in sunlight
beta = asind(sinB_all);               % beta angle [deg]
gamma = mean(beta);                   % mean cant angle for layout C [deg]

% Power obtained depending on solar array positioning strategies
P_A = P0 * flux;                              % ideal 2-axis tracking
P_B = P0 * flux .* max(0, sinB_all);          % body-fixed, normal = +h
P_C = P0 * flux .* max(0, cosd(beta - gamma));% 1 axis track about h, canted
P_D = P0 * flux .* max(0, cosd(beta));        % 1 axis track about h, no cant

% Energy per orbit
P_all = [P_A; P_B; P_C; P_D];               % 4xN [W]
strategy = ["A: 2-axis tracking", "B: body-fixed +h", ...
    "C: 1-axis canted",   "D: 1-axis no cant"];

% Power each strategy delivers while lit (no eclipse)[W]
P_lit = (P0 ./ dAU_all.^2) .* [ones(size(sinB_all));
                                max(0, sinB_all);
                                max(0, cosd(beta - gamma));
                                max(0, cosd(beta))];

tOrb  = eclStart;                           % orbit k starts at eclipse entry k
T_orb = diff(tOrb) / 60;                    % orbit durations [min]
T_lit = tOrb(2:end) - eclEnd(1:end-1);      % sunlit time of each orbit [s]
t_mid = (eclEnd(1:end-1) + tOrb(2:end)) / 2;% middle of each sunlit arc [s]
P_mid = interp1(T_all, P_lit.', t_mid).';   % lit power [W]
E_orb = P_mid .* T_lit / 3.6e6;             % energy per orbit [kWh]

% Results
fprintf("\nOrbits analysed: %d | mean duration %.2f min\n", numel(T_orb), mean(T_orb));
for k = 1:4
    fprintf("%-20s mean %.3f | min %.3f | max %.3f kWh/orbit\n", strategy(k), ...
        mean(E_orb(k,:)), min(E_orb(k,:)), max(E_orb(k,:)));
end

% Cross-check: total energy divided by number of orbits
E_tot = trapz(T_all, P_all, 2) / 3.6e6;                 % [kWh]
N_orb = (T_all(end) - T_all(1)) / (mean(T_orb)*60);     % orbits in the run
fprintf("Cross-check (total/N): %s kWh/orbit\n", mat2str(E_tot.'/N_orb, 4));

figure; plot(tOrb(1:end-1)/86400, E_orb, '.');
xlabel('Days since epoch'); ylabel('Energy per orbit [kWh]');
legend(strategy, 'Location', 'best'); grid on;

figure; plot(T_all/86400, beta);
xlabel('Days since epoch'); ylabel('\beta [deg]'); grid on;

%% Task 3(b)
clc; close all;

% Settings 
SelectOrbit = 2;                   
kep0 = elems(SelectOrbit,:);
a0   = kep0(1);                      

sim_days = 100;
stop_time  = start_time + days(sim_days);
sample_time  = 60;            
max_samples  = 20000;    
chunk_length = seconds(max_samples*sample_time);

n0    = sqrt(mu/a0^3);               % mean motion [rad/s]
T_orb = 2*pi/n0;                     % orbital period [s]
B     = Cd*area_drag/mass_sat;       % ballistic coefficient [m^2/kg]

% Chunk propagation with drag 
T_all = [];   % time since start_time [s]
a_all = [];   % osculating semi-major axis [m]
e_all = [];   % eccentricity 

kep = kep0;
t0  = start_time;
while t0 < stop_time
    t1 = min(t0 + chunk_length, stop_time);

    % Scenario for this chunk
    sc = satelliteScenario(t0, t1, sample_time);
    numericalPropagator(sc, ...
        'IncludeAtmosDrag', true, ...
        'GravitationalPotentialModel', "point-mass", ...
        'ODESolver', "ode45", ...
        'ODESet', odeset('RelTol',1e-8,'AbsTol',1e-8));

    sat = satellite(sc, kep(1), kep(2), kep(3), kep(4), kep(5), kep(6), ...
        "Name", "SAT_SEL", "OrbitPropagator", "numerical");
    physProps = physicalProperties(sat, 'Mass', mass_sat, ...
        'DragArea', area_drag, ...
        'DragCoefficient', Cd);

    % States of this chunk
    [r, v, t] = states(sat, "CoordinateFrame", "inertial");
    [a, ecc]  = ijk2keplerian(r, v);
    tsec = seconds(t - start_time);
    tsec = tsec(:);
    a = a(:);
    ecc = ecc(:);

    % The first sample of each new chunk repeats the last sample of the previous one
    if ~isempty(T_all)
        tsec = tsec(2:end);
        a = a(2:end);
        ecc = ecc(2:end);
    end
    T_all = [T_all; tsec];
    a_all = [a_all; a];
    e_all = [e_all; ecc];

    % Prepare next chunk propagation
    [r_end, v_end] = states(sat, t1, "CoordinateFrame", "inertial");
    [a, ecc, incl, RAAN, argp, nu] = ijk2keplerian(r_end, v_end);
    kep = [a ecc incl RAAN argp nu];

    fprintf("Done until %s\n", string(t1, 'dd-MMM-yyyy HH:mm'));
    t0 = t1;
    clear sc sat physProps r v t
end

% Remove short-period oscillations (orbit average)
window    = round(T_orb/sample_time);   % samples per orbit
a_mean = movmean(a_all, window);
idx = (window:numel(T_all)-window).';

% Secular decay rate and total loss in a 
p = polyfit(T_all(idx), a_mean(idx), 1);
adot = p(1);                       % [m/s]
T_win = sim_days*86400;            % [s]
Delta_a = adot*T_win;              % secular loss of a over 100 days [m]

fprintf('\n--- Propagation results (%d days) ---\n', sim_days);
fprintf('Decay rate of a : %.3f m/day\n', adot*86400);
fprintf('Total loss of a : %.2f m\n', Delta_a);
fprintf('Change in e : %.2e\n', e_all(end) - e_all(1));

% Station-keeping DeltaV
% Two tangential burns: dV = (n/2)*|Delta a|
dV_100   = 0.5*n0*abs(Delta_a);
% Equivalent continuous thrust to compensate drag
dVdt     = 0.5*n0*abs(adot);          % [m/s^2]
dV_orbit = dVdt*T_orb;
dV_day   = dVdt*86400;

fprintf('\n--- Station-keeping DeltaV ---\n');
fprintf('DeltaV per orbit    : %.3e m/s\n', dV_orbit);
fprintf('DeltaV per day      : %.3e m/s\n', dV_day);
fprintf('DeltaV in 100 days  : %.4f m/s\n', dV_100);

% Analytical cross-check 
% Density MATLAB's model is applying 
rho_eff = abs(adot)/(B*sqrt(mu*a0));
fprintf('\n--- Analytical check ---\n');
fprintf('Altitude            : %.1f km\n', (a0 - Re)/1e3);
fprintf('Effective density   : %.3e kg/m^3\n', rho_eff);

rho_tab = 6.65E-14; % Assignment guide density table value at 660 km (mean)
v0      = sqrt(mu/a0);
acc_drag  = 0.5*rho_tab*v0^2*B;               
dV_100_an = acc_drag*T_win;
da_rev    = -2*pi*rho_tab*a0^2*B;          % loss of a per revolution [m]
fprintf('Table density       : %.3e kg/m^3\n', rho_tab);
fprintf('Decay per orbit     : %.3f m\n', da_rev);
fprintf('DeltaV 100 days     : %.4f m/s\n', dV_100_an);

%  Operational window strategy
win_lim  = 100;            % station keeping window limit [m]
n_burns = ceil(abs(Delta_a)/win_lim);
dV_burn = 0.5*n0*win_lim;
t_burn  = win_lim/abs(adot)/86400;              % days between maneuvers
fprintf('\n--- Operational window (limit %.0f m) ---\n', win_lim);
fprintf('Maneuver every      : %.1f days\n', t_burn);
fprintf('Number of maneuvers : %d\n', n_burns);
fprintf('DeltaV per maneuver : %.4f m/s\n', dV_burn);

% Plots 
figure('Name', "Selected orbit - a decay");
plot(T_all/86400, (a_all - a0), 'Color', [0.7 0.7 1]); hold on
plot(T_all(idx)/86400, (a_mean(idx) - a0), 'b', 'LineWidth', 1.2)
plot(T_all(idx)/86400, (polyval(p, T_all(idx)) - a0), 'r--', 'LineWidth', 1.5)
xlabel('Time since epoch [days]'); ylabel('a - a_0 [m]'); grid on
legend('Osculating a', 'Orbit-averaged a', 'Linear fit', 'Location', 'southwest')
title('Selected orbit: semi-major axis decay due to drag (100 days)')

figure('Name', "Selected orbit - eccentricity");
plot(T_all/86400, e_all, 'b'); grid on
xlabel('Time since epoch [days]'); ylabel('e [-]')
title('Selected orbit: eccentricity with drag (100 days)')


%% Auxiliary functions


% Check whether an orbit is sun-synchronous given the initial orbital elements.
% If it is, also compute the Local Time of the Ascending Node (LTAN) at epoch.
function [isSSO, raanDot, sunDot, ltan, ltanApp] = isSunSynchronous(kep,...
                                            name, tol, epoch, Re, mu, J2)
n = sqrt(mu / kep(1)^3);
p = kep(1) * (1 - kep(2)^2);
K = 1.5 * n * J2 * (Re/p)^2;
raanDot = -K * cosd(kep(3)) * (180/pi) * 86400;   % secular RAAN rate [deg/day]
sunDot  = 360 / 365.2422;                     % Sun's mean motion [deg/day]
isSSO = abs(raanDot/sunDot - 1) <= tol;
ltanDrift = 4 * (raanDot - sunDot);   % LTAN drift [min/day] (1 deg = 4 min)
yesno = {'no', 'yes'};
fprintf('%s; dRAAN/dt = %+.5f deg/day vs Sun %.5f deg/day (%+.2f %%)\n', ...
    name, raanDot, sunDot, 100*(raanDot/sunDot - 1));
fprintf('%s; LTAN drift = %+.4f min/day (%+.1f min/year) -> sun-synchronous: %s (tol %.2g %%)\n', ...
    name, ltanDrift, 365.25*ltanDrift, yesno{isSSO + 1}, 100*tol);

ltan = NaN;
ltanApp = NaN;
if isSSO && nargin >= 7
    % Days since J2000.0 (2000-01-01 12:00 TT; the ~1 min UTC-TT offset is negligible here)
    if isdatetime(epoch)
        d = days(epoch - datetime(2000, 1, 1, 12, 0, 0));
    else
        d = epoch - 2451545.0;      % epoch given as a Julian Date
    end

    % Low-precision Sun position (Astronomical Almanac, ~0.01 deg accuracy)
    L   = mod(280.460 + 0.9856474*d, 360);      % Sun mean longitude [deg]
    g   = mod(357.528 + 0.9856003*d, 360);      % Sun mean anomaly [deg]
    lam = L + 1.915*sind(g) + 0.020*sind(2*g);  % Sun ecliptic longitude [deg]
    eps = 23.439 - 4e-7*d;                      % obliquity of the ecliptic [deg]
    alphaSun = mod(atan2d(cosd(eps)*sind(lam), cosd(lam)), 360);  % true Sun RA [deg]

    % LTAN: the Sun is at local noon (12 h) where the node's RA equals the Sun's RA,
    % and every 15 deg of RA east of the Sun adds 1 h of local time.
    ltan    = mod(12 + (kep(4) - L)/15, 24);          % mean Sun  (RA of mean Sun = L)
    ltanApp = mod(12 + (kep(4) - alphaSun)/15, 24);   % true Sun

    fprintf('%s | Sun RA at epoch = %.3f deg (mean Sun %.3f deg)\n', name, alphaSun, L);
    fprintf('%s | LTAN at epoch = %s mean solar time (%s apparent solar time)\n', ...
        name, hhmm(ltan), hhmm(ltanApp));
end
end

% Format decimal hours as hh:mm
function s = hhmm(h)
m = mod(round(h*60), 1440);
s = sprintf('%02d:%02d', floor(m/60), mod(m, 60));
end

% Enter NaN between jumps of more than 180 deg
function [latP, lonP] = breakAtDateline(lat, lon)
    k = find(abs(diff(lon)) > 180);       
    lonP = lon;  latP = lat;
    for j = flip(k)                            
        lonP = [lonP(1:j), NaN, lonP(j+1:end)];
        latP = [latP(1:j), NaN, latP(j+1:end)];
    end
end

% Appends the intervals of a new chunk to the global list. If the previous
% chunk ended while in contact and the new chunk starts in contact, it is
% the same pass split by the boundary, so the two pieces are joined.
function [allS, allE] = appendIntervals(allS, allE, newS, newE, chunkStart, tol)
    newS = newS(:);  newE = newE(:);
    if ~isempty(allE) && ~isempty(newS) && ...
            abs(allE(end) - chunkStart) < tol && ...   % previous pass reached the boundary
            abs(newS(1)   - chunkStart) < tol          % new pass starts at the boundary
        allE(end) = newE(1);      % extend the previous pass
        newS(1) = [];  newE(1) = [];
    end
    allS = [allS; newS];
    allE = [allE; newE];
end

function [ha, hp, T] = getOrbitCharacteristics(kep, Re, mu)
    ha = (kep(1)*(1+kep(2)) - Re)/1000; % Initial apogee altitude [km]
    hp = (kep(1)*(1-kep(2)) - Re)/1000; % Initial perigee altitude [km]
    n0    = sqrt(mu/kep(1)^3);               
    T = 2*pi/n0;                        % Initial orbital period
end
