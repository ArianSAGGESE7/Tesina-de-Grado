


% close all;clear all
% % datos = loadAlmanacGps("almanac.yuma.week0100.061440.txt");
% datos = loadAlmanacGps("usat.yuma.txt");
% 
% % Definimos satélite a utilizar
% 
% Sv=1;
% 
% datosPRN=datos(Sv); % Datos del satélite Sv
% 
% 
% 
% [pos, vel] = Orbital2ECEF(t, datosPRN.toe, ...
%     datosPRN.sqa^2, 0, 0, ...
%     3.9860e+14, 0, 0, ...
%     datosPRN.m0, datosPRN.e, datosPRN.om0, ...
%     datosPRN.i0, 0, ...
%     0, 0, 0, 0, 0, 0, ...
%     7.2921151467e-5, datosPRN.omd, datosPRN.w,0,datosPRN.af0,datosPRN.af1);



function [pos, vel] = Orbital2ECEF(t, toe, ...
    ARef, deltaA, ADot, ...
    mu, deltan0, deltan0Dot, ...
    M0, e0, omega0, ...
    i0, iDot, ...
    Cis, Cic, Crs, Crc, Cus, Cuc, ...
    OmegaEDot, OmegaRefDot, Omega0,deltaOmegaDot,af0,af1)

% Función de parámetros orbitales a coordenadas ECEF


% ORBITPARAMETERSTOECEF

%   Convert orbital parameters to satellite positions
%   and velocities in Earth-Centered Earth-Fixed (ECEF)
%
%   All inputs should be of the same size S-by-1-by-T, where S is the
%   number of satellites and T is the number of time samples.
%
%   All outputs are of size S-by-1-by-T, where S is the number of
%   satellites and T is the number of time samples.
%
%   gpsWeek       - Current GPS week number
%   t             - Current time in GPS week
%
%   Ephemeris parameters
%
%   weekNum       - Week number
%   toe           - Reference time of week (s)
%   ARef          - Reference semi-major axis (m) (Constant)
%   deltaA        - Semi-major axis difference at reference time (m)
%   ADot          - Change rate in semi-major axis (m/s)
%   mu            - WGS 84 value of Earth's gravitational constant
%                   (m^3/s^2) (Constant)
%   deltan0       - Mean motion difference from computed value at reference
%                   time (rad/s)
%   deltan0Dot    - Rate of mean motion difference from computed value
%                   (rad/s^2)
%   M0            - Mean anomaly at reference time (rad)
%   e0            - Eccentricity
%   omega0        - Argument of perigee (rad)
%   i0            - Inclination angle at reference time (rad)
%   iDot          - Rate of inclination angle (rad/s)
%   Cis           - Amplitude of the sine harmonic correction term to the
%                   angle of inclination (rad)
%   Cic           - Amplitude of the cosine harmonic correction term to the
%                   orbit radius (m)
%   Crs           - Amplitude of the sine harmonic correction term to the
%                   orbit radius (m)
%   Crc           - Amplitude of the cosine harmonic correction term to the
%                   orbit radius (m)
%   Cus           - Amplitude of the sine harmonic correction term to the
%                   argument of latitude (rad)
%   Cuc           - Amplitude of the cosine harmonic correction term to the
%                   argument of latitude (rad)
%   OmegaEDot     - WGS 84 value of Earth's rotation rate (rad/s)
%                   (Constant)
%   OmegaRefDot   - Reference rate of right of ascension (rad/s) (Constant)
%   Omega0        - Longitude of ascending node of orbit plane at weekly
%                   epoch (rad)
%   deltaOmegaDot - Rate of right ascension difference (rad/s)



% Time from ephemeris reference time
tSys = t;
tk =  tSys - toe - (af0 + af1*(tSys-toe)); % Correccion de fenomenos de clk

% Semi-Major Axis at reference time
A0 = ARef + deltaA;
% Semi-Major Axis
Ak = A0 + ADot.*tk;

% Computed Mean Motion
n0 = sqrt(mu ./ A0.^3);
% Mean motion difference from computed value
deltanA = deltan0 + 0.5.*deltan0Dot.*tk;
% Corrected Mean Motion
nA = n0 + deltanA;

% Mean Anomaly
Mk = M0 + nA.*tk;
% Eccentric Anomaly
Ek = solveKeplerEquationEccentricAnomaly(e0, Mk);
cosEk = cos(Ek);
sinEk = sin(Ek);
% True Anomaly
vk = atan2((sqrt(1 - e0.^2) .* sinEk) ./ (1 - e0 .* cosEk), ...
    (cosEk - e0) ./ (1 - e0 .* cosEk));
% Argument of Latitude
Phik = vk + omega0;
sin2Phik = sin(2*Phik);
cos2Phik = cos(2*Phik);

% Calculate second harmonic perturbations
% Argument of Latitude Correction
deltauk = Cus .* sin2Phik + Cuc .* cos2Phik;
% Radial Correction
deltark = Crs .* sin2Phik + Crc .* cos2Phik;
% Inclination Correction
deltaik = Cis .* sin2Phik + Cic .* cos2Phik;

% Corrected Argument of Latitude
uk = Phik + deltauk;
cosuk = cos(uk);
sinuk = sin(uk);
% Corrected Radius
rk = Ak .* (1 - e0 .* cosEk) + deltark;
% Corrected Inclination
ik = i0 + iDot .* tk + deltaik;
cosik = cos(ik);
sinik = sin(ik);

% Positions in orbital plane
xkPrime = rk .* cosuk;
ykPrime = rk .* sinuk;

% Rate of Right Ascension
OmegaDot = OmegaRefDot + deltaOmegaDot;
% Corrected Longitude of Ascending Node
Omegak = Omega0 + (OmegaDot - OmegaEDot) .* tk - OmegaEDot .* toe;
cosOmegak = cos(Omegak);
sinOmegak = sin(Omegak);

% ECEF Positions of Satellites
xk = xkPrime .* cosOmegak - ykPrime .* cosik .* sinOmegak;
yk = xkPrime .* sinOmegak + ykPrime .* cosik .* cosOmegak;
zk = ykPrime .* sinik;

% Concatenate outputs. S-by-3-by-T output matrix where S is the number of
% satellites and T is the number of time samples.
pos = [xk, yk, zk];

% Eccentric Anomaly derivative
EkDot = nA ./ (1 - e0 .* cosEk);
% Argument of Latitude derivative
PhikDot = (sin(vk) ./ sinEk) .* EkDot;

% Argument of Latitude Correction derivative
deltaukDot = 2 * (Cus .* cos2Phik - Cuc .* sin2Phik) .* PhikDot;
% Radial Correction derivative
deltarkDot = 2 * (Crs .* cos2Phik - Crc .* sin2Phik) .* PhikDot;
% Inclination Correction derivative
deltaikDot = 2 * (Cis .* cos2Phik - Cic .* sin2Phik) .* PhikDot;

% Corrected Radius derivative
rkDot = (Ak .* e0 .* sinEk) .* EkDot + deltarkDot;
% Corrected Argument of Latitude derivative
ukDot = PhikDot + deltaukDot;
% Corrected Inclination derivative
ikDot = iDot + deltaikDot;

% Velocities in orbital plane
xkPrimeDot = rkDot .* cosuk - rk .* ukDot .* sinuk;
ykPrimeDot = rkDot .* sinuk + rk .* ukDot .* cosuk;

% Corrected Longitude of Ascending Node derivative
OmegakDot = OmegaDot - OmegaEDot;

% ECEF Velocities of Satellites
xkDot = xkPrimeDot .* cosOmegak - ykPrimeDot .* cosik .* sinOmegak ...
    + ikDot .* ykPrime .* sinik .* sinOmegak - OmegakDot .* yk;
ykDot = xkPrimeDot .* sinOmegak + ykPrimeDot .* cosik .* cosOmegak ...
    - ikDot .* ykPrime .* sinik .* cosOmegak + OmegakDot .* xk;
zkDot = ykPrimeDot .* sinik + ikDot .* ykPrime .* cosik;

% Concatenate outputs. S-by-3-by-T output matrix where S is the number of
% satellites and T is the number of time samples.
vel = [xkDot, ykDot, zkDot];
end

function Ek = solveKeplerEquationEccentricAnomaly(e0, Mk)
% Solve Kepler's equation for Eccentric Anomaly
% Mk = Ek - e0*sin(Ek)
%
% Inputs:
% e0 Eccentricity
% Mk Mean Anomaly (rad)
%
% Outputs:
% Ek Eccentric Anomaly (rad)

numIters = 22;  % Empirical number to achieve millimeter accuracy.

Ek = Mk + (e0 .* sin(Mk)) ./ (1 - sin(Mk + e0) + sin(Mk));
for iter = 1:numIters
    Ek = Mk + e0 .* sin(Ek);
end

end


