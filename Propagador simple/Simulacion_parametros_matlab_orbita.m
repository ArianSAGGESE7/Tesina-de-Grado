%Simulación Matlab 


startTime = datetime(2020,5,10);
stopTime = startTime + days(5);
sampleTime = 60;                                       % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);

earthAngularVelocity = 0.0000729211585530;                                             % rad/s
orbitalPeriod = 2*pi/earthAngularVelocity;                                             % seconds
earthStandardGravitationalParameter = 398600.4418e9;                                   % m^3/s^2
semiMajorAxis = (earthStandardGravitationalParameter*((orbitalPeriod/(2*pi))^2))^(1/3);

eccentricity = 0;
inclination = 60;                  % degrees
rightAscensionOfAscendingNode = 0; % degrees
argumentOfPeriapsis = 0;           % degrees
trueAnomaly = 0;                   % degrees

sat = satellite(sc,semiMajorAxis,eccentricity,inclination,rightAscensionOfAscendingNode,...
        argumentOfPeriapsis,trueAnomaly,"OrbitPropagator","two-body-keplerian","Name","GEO Sat");
v = satelliteScenarioViewer(sc);

%% Ground track

leadTime = 2*24*3600;                                          % seconds
trailTime = leadTime;
gt = groundTrack(sat,"LeadTime",leadTime,"TrailTime",trailTime);


play(sc);