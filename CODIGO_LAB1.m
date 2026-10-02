%%%%LABORATORY SESSION 1%%%%
%% SAT A

clc;clear;close all;

% %COE
a_a = 26556e3; %semimayor axis
e_a = 0.70399; %eccenntricity
i_a = 63.4; %inclination
raan_a = 90; %right ascension of the ascending node
aop_a = 270; %argument of periapsis
ta_a = 0; %true anomaly

%Select start date, end date and sample time 
start_time = datetime(2013,11,10,23,03,03.95);
stop_time = start_time+days(365); %propagate for 10 day
sample_time=600; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

satA = satellite(sc,a_a,e_a,i_a,raan_a,aop_a,ta_a,"Name","SAT_A", ...
    "OrbitPropagator","two-body-keplerian");

%Add a sensor 
g = gimbal(satA); %CHANGE HERE WHICH SAT YOU WANT TO STUDY
pointAt(g,'nadir');
camera= conicalSensor(g,MaxViewAngle=20); %Add a camera with gimbal g and a fov of 20 deg
fieldOfView(camera);
groundTrack(satA,LeadTime=86400*30)

%Add ground station
lat = 40;
lon = -3.7033;
min_elev = 20; %minimum elevation angle 
gs = groundStation(sc,"Latitude",lat,"Longitude",lon,"MaskElevationAngle",min_elev,"Name","GS");

% satelliteScenarioViewer(sc,"Dimension","3D");
% satelliteScenarioViewer(sc,"Dimension","2D");

ac = access(satA,gs);
intvls = accessIntervals(ac)

window = zeros(length(intvls.StartTime));

for i = 1:length(window)
    window(i) = minutes(intvls.EndTime(i) - intvls.StartTime(i));
end

avgWindow = mean(window);

fprintf("Average window duration= %.2f min\n",avgWindow(1));

fprintf("Number of contacts in a year= %.2f\n", length(intvls.IntervalNumber))

fprintf("Frequency of contact= %.2f/day\n",length(intvls.IntervalNumber)/365);


%% SAT B

clc;clear;close all;

%COE
a_b = 7039.5e3; %semimayor axis
e_b = 0.0001; %eccenntricity
i_b = 98.074; %inclination
raan_b = 287; %right ascension of the ascending node
aop_b = 260.31; %argument of periapsis
ta_b = 0; %true anomaly

%Select start date, end date and sample time 
start_time = datetime(2013,11,10,23,03,03.95);
stop_time = start_time+days(220); %propagate for 10 day
sample_time=300; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

satB = satellite(sc,a_b,e_b,i_b,raan_b,aop_b,ta_b,"Name","SAT_B", ...
    "OrbitPropagator","SGP4");

%Add a sensor 
g = gimbal(satB); %CHANGE HERE WHICH SAT YOU WANT TO STUDY
pointAt(g,'nadir');
camera= conicalSensor(g,MaxViewAngle=20); %Add a camera with gimbal g and a fov of 20 deg
fieldOfView(camera);
% groundTrack(satB,LeadTime=86400*5)

%Add ground station
lat = 40;
lon = -3.7033;
min_elev = 20; %minimum elevation angle 
gs = groundStation(sc,"Latitude",lat,"Longitude",lon,"MaskElevationAngle",min_elev);

% satelliteScenarioViewer(sc,"Dimension","3D");
% satelliteScenarioViewer(sc,"Dimension","2D");

ac = access(satB,gs);
intvls = accessIntervals(ac);

window = zeros(length(intvls.StartTime));

for i = 1:length(window)
    window(i) = minutes(intvls.EndTime(i) - intvls.StartTime(i));
end

avgWindow = mean(window);

fprintf("Average window duration= %.2f min\n",avgWindow(1));

fprintf("Number of contacts in a year= %.2f\n", length(intvls.IntervalNumber))

fprintf("Frequency of contact= %.2f/day\n",length(intvls.IntervalNumber)/365);


%% SAT C

clc;clear;close all;

%COE
a_c = 7202.5e3; %semimayor axis
e_c = 0.02; %eccenntricity
i_c = 42.2; %inclination
raan_c = 121; %right ascension of the ascending node
aop_c = 15; %argument of periapsis
ta_c = 0; %true anomaly

%Select start date, end date and sample time 
start_time = datetime(2013,11,10,23,03,03.95);
stop_time = start_time+days(365); %propagate for 10 day
sample_time=600; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);

satC = satellite(sc,a_c,e_c,i_c,raan_c,aop_c,ta_c,"Name","SAT_C", ...
    "OrbitPropagator","two-body-keplerian");

%Add a sensor 
g = gimbal(satC); %CHANGE HERE WHICH SAT YOU WANT TO STUDY
pointAt(g,'nadir');
camera= conicalSensor(g,MaxViewAngle=20); %Add a camera with gimbal g and a fov of 20 deg
fieldOfView(camera);
groundTrack(satC,LeadTime=86400*10)

%Add ground station
lat = 40;
lon = -3.7033;
min_elev = 20; %minimum elevation angle 
gs = groundStation(sc,"Latitude",lat,"Longitude",lon,"MaskElevationAngle",min_elev);

% satelliteScenarioViewer(sc,"Dimension","3D");
% satelliteScenarioViewer(sc,"Dimension","2D");

ac = access(satC,gs);
intvls = accessIntervals(ac)

window = zeros(length(intvls.StartTime));

for i = 1:length(window)
    window(i) = minutes(intvls.EndTime(i) - intvls.StartTime(i));
end

avgWindow = mean(window);

fprintf("Average window duration= %.2f min\n",avgWindow(1));

fprintf("Number of contacts in a year= %.2f\n", length(intvls.IntervalNumber))

fprintf("Frequency of contact= %.2f/day\n",length(intvls.IntervalNumber)/365);


%% SAT D

%COE
a_d = 7202.5e3; %semimayor axis
e_d = 0.0009826; %eccenntricity
i_d = 96.5717; %inclination
raan_d = 344.5256; %right ascension of the ascending node
aop_d = 296.2811; %argument of periapsis
ta_d = 0; %true anomaly

%Select start date, end date and sample time 
start_time = datetime(2013,11,10,23,03,03.95);
stop_time = start_time+days(360); %propagate for 10 day
sample_time=60; %Sample every 60 seconds
%Create Satellite scenario
sc = satelliteScenario(start_time,stop_time,sample_time);
satD = satellite(sc,a_d,e_d,i_d,raan_d,aop_d,ta_d,"Name","SAT_D", ...
    "OrbitPropagator","two-body-keplerian");

%Add a sensor 
g = gimbal(satC); %CHANGE HERE WHICH SAT YOU WANT TO STUDY
pointAt(g,'nadir');
camera= conicalSensor(g,MaxViewAngle=20); %Add a camera with gimbal g and a fov of 20 deg
fieldOfView(camera);
groundTrack(satC,LeadTime=86400*5)

%Add ground station
lat = 40;
lon = -3.7033;
min_elev = 20; %minimum elevation angle 
gs = groundStation(sc,"Latitude",lat,"Longitude",lon,"MaskElevationAngle",min_elev);

% satelliteScenarioViewer(sc,"Dimension","3D");
% satelliteScenarioViewer(sc,"Dimension","2D");

ac = access(satD,gs);
intvls = accessIntervals(ac)

window = zeros(length(intvls.StartTime));

for i = 1:length(window)
    window(i) = minutes(intvls.EndTime(i) - intvls.StartTime(i));
end

avgWindow = mean(window)

fprintf("Average window duration= %.2f min\n",avgWindow(1));

fprintf("Number of contacts in a year= %.2f\n", length(intvls.IntervalNumber))

fprintf("Frequency of contact= %.2f/day\n",length(intvls.IntervalNumber)/365);
