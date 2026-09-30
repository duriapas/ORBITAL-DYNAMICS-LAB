%% Creating a satellite scenario 
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 10 day
sample_time=200 ; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

% Add a satellite to the satellite scenario 
%We can specify orbital elements or TLE file, which encodes the orbital
%elements for objects orbiting around earth 
%https://es.mathworks.com/help/aerotbx/ug/satelliteScenario-key-concepts.html#mw_e08739b4-c5da-4983-898b-56e18cc71f87
%There are elements of tracked objects (real satellite, this is extremely 
%helpful)
%lte files available at https://www.space-track.org/
tle_file = "../TLE/GALILEO5.tle";
sat1 = satellite(sc,tle_file,"Name","satTwoBodyKeplerian","OrbitPropagator","two-body-keplerian");

% Add a ground track for satellite 1 
groundTrack(sat1,LeadTime=86400*5)

%View satellite 
% satelliteScenarioViewer(sc,'Dimension','3D','CameraReferenceFrame','Inertial');
sat1.MarkerColor = "blue";
sat1.Orbit.LineColor = "green";
satelliteScenarioViewer(sc,'Dimension','2D');





%% Satellite defined by orbital elements
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 10 day
sample_time=200 ; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
a = 2.66e7; %semimajor axis in m
e = 0.74; %Eccentricity
i = 63.43; %Inclination
raan = 0; %Right ascension of the ascending node 
aop = 270; %argument of periapsis 
ta = 0; %True anomaly 


sat2 = satellite(sc,a,e,i,raan,aop,ta,"Name","sat-sgp4", ...
    "OrbitPropagator","sgp4");

% Add a ground track for satellite 1 
groundTrack(sat2,LeadTime=86400*5)
sat2.MarkerColor = "red";
sat2.Orbit.LineColor = "white";
%View satellite 
v = satelliteScenarioViewer(sc,'Dimension','2D');



%% Groundtrack of GEO satellite 
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 10 day
sample_time=200 ; %Sample every 200 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
a = 42164e3 ; %semimajor axis in m
e = 0.0; %Eccentricity
i = 30 ; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 0; %argument of periapsis 
ta = 0; %True anomaly 


sat2 = satellite(sc,a,e,i,raan,aop,ta,"Name","sat-sgp4", ...
    "OrbitPropagator","two-body-keplerian");

% Add a ground track for satellite 1 
groundTrack(sat2,LeadTime=86400*5)

%View satellite 
sat2.MarkerColor = "red";
sat2.Orbit.LineColor = "white";
satelliteScenarioViewer(sc,'Dimension','3D');


%% Groundtrack of GEO satellite 
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 10 day
sample_time=200 ; %Sample every 200 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
a = 42164e3 ; %semimajor axis in m
e = 0.0; %Eccentricity
i = 30 ; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 0; %argument of periapsis 
ta = 0; %True anomaly 

     numericalPropagator(sc,...
    'IncludeAtmosDrag', false, ...
    'GravitationalPotentialModel', "spherical-harmonics", ...
    'SphericalHarmonicModel','EGM2008',...
    'ODESolver', "ode45", ... 
    'ODESet', odeset('RelTol', 1e-8, 'AbsTol', 1e-8)); 


sat2 = satellite(sc,a,e,i,raan,aop,ta,"Name","sat-sgp4", ...
    "OrbitPropagator","numerical");

% Add a ground track for satellite 1 
groundTrack(sat2,LeadTime=86400*5)
sat2.MarkerColor = "red";
sat2.Orbit.LineColor = "white";
%View satellite 
v = satelliteScenarioViewer(sc,'Dimension','3D');


%% Satellite with drag, propagate orbit
clc;clear;
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(3.8); %propagate for 1 day

%Create Satellite scenario
a = 6533.137e3; %semimajor axis in m
e = 0.0; %Eccentricity
i = 0; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 0; %argument of periapsis 
ta = 0; %True anomaly 

%Satellite properties 
mass_sat  = 9e3;      % Satellite mass [kg]
area_drag = 1.1;        % Cross-sectional area S [m^2]
Cd        = 3.0;       % Drag coefficient [-]


% % Define spacecraft physical properties
physProps = Aero.spacecraft.PhysicalProperties ('Mass', mass_sat, ...
    'DragArea', area_drag, ...
    'DragCoefficient', Cd);
numOpts = Aero.spacecraft.NumericalPropagatorOptions(...
    'IncludeAtmosDrag', true, ...
    'GravitationalPotentialModel', "spherical-harmonics", ...
    'SphericalHarmonicModel','EGM96',...
    'ODESolver', "ode45", ... 
    'ODESet', odeset('RelTol', 1e-8, 'AbsTol', 1e-8)); 
% Define time step and time vector 
prop_sample_time = 10; % [s]
num_steps = round(seconds(stop_time - start_time) / prop_sample_time);
time_vec  = start_time + seconds((0:num_steps) * prop_sample_time);


[pos_icrf, vel_icrf] = propagateOrbit(time_vec, a, e, i, raan, aop, ta, ...
    "PropModel", "numerical", ...
    "NumericalPropagatorOptions", numOpts, ...
    "PhysicalProperties", physProps);


r=zeros(1,length(time_vec));
for i=1:length(time_vec)
    r(i) = (sum(pos_icrf(:,i).^2))^0.5;
end 
plot(time_vec,(r-6371e3)/1000)

%% Satellite with drag, satellitescenario
clc;clear;
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(0.60); %propagate for 1 day
sample_time = 2;
sc = satelliteScenario(start_time,stop_time,sample_time,AutoSimulate=false);

%Create Satellite scenario
a = 6633.137e3; %semimajor axis in m
e = 0.0; %Eccentricity
i = 0; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 270; %argument of periapsis 
ta = 0; %True anomaly 

%Satellite properties 
mass_sat  = 2e3;      % Satellite mass [kg]
area_drag = 15;        % Cross-sectional area S [m^2]
Cd        = 3.0;       % Drag coefficient [-]


% % Define spacecraft physical properties

numOpts = numericalPropagator(sc,...
    'IncludeAtmosDrag', true, ...
    'GravitationalPotentialModel', "spherical-harmonics", ...
    'SphericalHarmonicModel','EGM2008',...
    'ODESolver', "ode45", ... 
    'ODESet', odeset('RelTol', 1e-8, 'AbsTol', 1e-8)); 
% Define time step and time vector 
sat = satellite(sc,a,e,i,raan,aop,ta,"Name","sat", ...
    "OrbitPropagator","numerical");

physProps = physicalProperties (sat,'Mass', mass_sat, ...
    'DragArea', area_drag, ...
    'DragCoefficient', Cd);

[pos_icrf,vel_icrf,time_vec]=states(sat);
% [a1,~,~,~,~,~,~,~,~] = ijk2keplerian(pos_icrf, vel_icrf);
for i=1:length(time_vec)
    r(i) = (sum(pos_icrf(:,i).^2))^0.5;
end 
plot(time_vec,(r-6371e3)/1000)




%% Satellite with camera and FOV
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 1 day
sample_time=200 ; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
a = 2.66e7; %semimajor axis in m
e = 0.74; %Eccentricity
i = 63.43; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 270; %argument of periapsis 
ta = 0; %True anomaly 


sat2 = satellite(sc,a,e,i,raan,aop,ta,"Name","sat-sgp4", ...
    "OrbitPropagator","sgp4");

%Add a sensor 
g = gimbal(sat2);
pointAt(g,'nadir');
camera= conicalSensor(g,MaxViewAngle=20); %Add a camera with gimbal g and a fov of 20 deg



% Get ground station contact with the field of view of the camera
%View satellite and FOV
fieldOfView(camera);

% Add a ground track for satellite 1 
groundTrack(sat2,LeadTime=86400*5)

%View satellite 
v = satelliteScenarioViewer(sc,'Dimension','3D');





%% Satellite with ground station 
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 1 day
sample_time=200 ; %Sample every 60 seconds
%Creat Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
%Ground station is defined by latitude and longitude altitude and minimum elevation
%angle 

lat = 10;
lon = -30;
min_elev = 20; %minimum elevation angle 
gs = groundStation(sc,"Latitude",lat,"Longitude",lon,"MaskElevationAngle",min_elev);
% Add satellite 
tle_file = "../TLE/GALILEO5.tle";
sat1 = satellite(sc,tle_file,"Name","satTwoBodyKeplerian", ...
    "OrbitPropagator","sgp4");
sat1.MarkerColor = "blue";
sat1.Orbit.LineColor = "green";

% Add a ground track for satellite 1 
groundTrack(sat1,LeadTime=86400*10)

% Get ground station contact information 
ac = access(sat1,gs);
intvls = accessIntervals(ac)
%View satellite 
play(sc);


%% Photograph object and calculate revisit time 

%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(10); %propagate for 1 day
sample_time=200 ; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

%Location of Madrid 

lat = 40.4168;
lon = -3.7038;
altitude = 600;
gs = groundStation(sc,"Latitude",lat,"Longitude",lon,"Altitude",altitude);

%Satellite defined by orbital elements
a = 1e7; %semimajor axis in m
e = 0.0; %Eccentricity
i = 30; %Inclination
raan = 270; %Right ascension of the ascneding node 
aop = 0; %argument of periapsis 
ta = 0; %True anomaly 


sat1 = satellite(sc,a,e,i,raan,aop,ta,"Name","satTwoBodyKeplerian", ...
    "OrbitPropagator","sgp4");
sat1.MarkerColor = "blue";
sat1.Orbit.LineColor = "green";

% Add gimbal for control of sensor 
g = gimbal(sat1);
%The gimbal can be steered indipendently of the satellite position 

%Track the location to be photographed with the gimbal 
pointAt(g,gs)
%We can also point at 'nadir' or give attitude tables

%Add a sensor 
camera= conicalSensor(g,MaxViewAngle=20); %Add a camera with gimbal g and a fov of 20 deg



% Get ground station contact with the field of view of the camera
ac = access(camera,gs);
%Check contacts 
intvls = accessIntervals(ac);
%View satellite and FOV
fieldOfView(camera);
play(sc);

%Compute max revisit time 
t_start = intvls.StartTime;
t_end = intvls.EndTime;
revisitTimes = hours(t_start(2:end) - t_end(1:end-1));
maxRevisitTime = max(revisitTimes); %In hours

%% Eclipse 
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(5); %propagate for 1 day
sample_time=20 ; %Sample every 60 seconds
%Create Satellite scenario

sc = satelliteScenario(start_time,stop_time,sample_time);
a = 1e7; %semimajor axis in m
e = 0.0; %Eccentricity
i = 20; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 270; %argument of periapsis 
ta = 0; %True anomaly 


sat1 = satellite(sc,a,e,i,raan,aop,ta,"Name","eclipse", ...
    "OrbitPropagator","sgp4");
%Compute eclipse
eclSat = eclipse(sat1);
v = satelliteScenarioViewer(sc,'Dimension','2D');
intvls = eclipseIntervals(eclSat)

%% Effect of J2 
%Select start date, end date and sample time 
start_time = datetime(2026,7,8,12,0,0);
stop_time = start_time+days(100); %propagate for 1 day
sample_time=200 ; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
a = 2.66e7; %semimajor axis in m
e = 0.74; %Eccentricity
i = 63.43; %Inclination
raan = 0; %Right ascension of the ascneding node 
aop = 270; %argument of periapsis 
ta = 0; %True anomaly 

sat1 = satellite(sc,a,e,i,raan,aop,ta,"Name","satTwoBodyKeplerian", ...
    "OrbitPropagator","sgp4");

sat2 = satellite(sc,a,e,30,raan,aop,ta,"Name","satTwoBodyKeplerian", ...
    "OrbitPropagator","sgp4");

g1=groundTrack(sat2,LeadTime=86400*5);
g2=groundTrack(sat1,LeadTime=86400*5);
g1.LeadLineColor = "red"
g2.LeadLineColor = "blue"

%View satellite 
sat2.MarkerColor = "red";
sat2.Orbit.LineColor = "white";
v = satelliteScenarioViewer(sc,'Dimension','2D');


% Extract evolution of orbital elements
[pos_1,vel_1]   = states(sat1);
[pos_2,vel_2] = states(sat2);
[~,~,~,~,argp_1,~,~,~,~] = ijk2keplerian(pos_1, vel_1);
[~,~,~,~,argp_2,~,~,~,~] = ijk2keplerian(pos_2, vel_2);
plot(argp_1);
hold on
plot(argp_2);