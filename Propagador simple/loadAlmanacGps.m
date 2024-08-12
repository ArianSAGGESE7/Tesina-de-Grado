function al = loadAlmanacGps(NameOfFile)
% Esta funci´on se utiliza para levantar los datos desde un archivo de
% almanaque YUMA
% Estructura donde se almacenan los par´ametros.
ALM = struct(...
'ID',... % PRN of the SVN
0,... %
'Health',... % 000=usable
0,... %
'e',... % Eccentricity: This shows the amount of the orbit
0,... % deviation from circular (orbit). It is the distance
... % between the foci divided by the length of the
... % semi-major axis (our orbits are very circular).
'toe',... % Time of Applicability: The number of seconds in the
0,... % orbit when the almanac was generated. [s]
... % Kind of a time tag.
'i0',... % Orbital Inclination: The angle to which the SV orbit
0,... % meets the equator (GPS is at approx. 55 degrees).
... % Roughly, the SV's orbit will not rise above approx.
... % 55 degrees latitude. The number is part of an
... % equation: #= pi/180 = the true inclination. [rad]
'omd',... % Rate of Right Ascension: Rate of change in the
0,... % measurement of the angle of right ascension as
... % defined in the Right Ascension mnemonic. [r/s]
'sqa',... % SQRT(A): Square Root of Semi-Major Axis: This is
0,... % defined as the measurement from the center of the
... % orbit to either the point of apogee or the point of
... % perigee. [m 1/2]
'om0',... % Right Ascension at Time of Almanac (TOA): Geographic
0,... % longitude of the ascending node of the orbit plane at
... % the weekly epoch. [rad]
'w',... % Argument of Perigee: An angular measurement along
0,... % the orbital path measured from the ascending node to
... % the point of perigee, measured in the direction of
... % the SV's motion. [rad]
'm0',... % Mean Anomaly: Angle (arc) traveled past the
0,... % longitude of ascending node (value= 0-180 degrees or
... % 0-negative 180 degrees). If the value exceeds 180
... % degrees, subtract 360 degrees to find the mean
... % anomaly. When the SV has passed perigee and heading
... % towards apogee, the mean anomaly is positive. After
... % the point of apogee, the mean anomaly value will be
... % negative to the point of perigee. [rad]
'af0',... % Af(0): SV clock bias in seconds. [s]
0,... %
'af1',... % Af(1): SV clock Drift in seconds per seconds [s/s]
0,... %
'week',... % GPS week (0000-1024), every 7 days since 22 Aug 1999
0 ... %
);
al(1:32) = ALM;
fp = fopen(NameOfFile,'r');
% N´umero de palabras que forman la definici´on de cada uno de los campos a leer
nOfWords = [1,1,1,3,2,4,3,4,3,2,1,1,1];
% Elimino la l´?nea de texto
while(~isempty(fscanf(fp,' %s',7)))
for k = 1:13
textString = fscanf(fp,' %s',nOfWords(k));
switch(textString)
case 'ID:'
svid = str2num(fscanf(fp,' %s',1));
al(svid).ID = svid;
case 'Health:'
al(svid).Health = str2num(fscanf(fp,' %s',1));
case 'Eccentricity:'
al(svid).e = str2num(fscanf(fp,' %s',1));
case 'TimeofApplicability(s):'
al(svid).toe = str2num(fscanf(fp,' %s',1));
case 'OrbitalInclination(rad):'
al(svid).i0 = str2num(fscanf(fp,' %s',1));
case 'RateofRightAscen(r/s):'
al(svid).omd = str2num(fscanf(fp,' %s',1));
case 'SQRT(A)(m1/2):'
al(svid).sqa = str2num(fscanf(fp,' %s',1));
case 'RightAscenatWeek(rad):'
al(svid).om0 = str2num(fscanf(fp,' %s',1));
case 'ArgumentofPerigee(rad):'
al(svid).w = str2num(fscanf(fp,' %s',1));
case 'MeanAnom(rad):'
al(svid).m0 = str2num(fscanf(fp,' %s',1));
case 'Af0(s):'
al(svid).af0 = str2num(fscanf(fp,' %s',1));
case 'Af1(s/s):'
al(svid).af1 = str2num(fscanf(fp,' %s',1));
case 'week:'
al(svid).week = str2num(fscanf(fp,' %s',1));
otherwise
disp('Error');
disp(textString);
end;
end;
end
fclose(fp);