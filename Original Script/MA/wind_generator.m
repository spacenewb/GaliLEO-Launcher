function [uw,vw] = wind_generator(z,Y,site)
% 
uw =0;
vw =0;

% 
% if z < 0
%     z = 0;
% end
% 
% if z < 70000
%     % diifferentiate for para simulation
%     xnorth = cosd(site.azimuth).*Y(1);
%     yeast = sind(site.azimuth).*Y(1);
%     [lat,lon,h] = ned2geodetic(xnorth, yeast, -z, site.lat0, site.lon0, -site.z0, wgs84Ellipsoid);
%     
%     if h > 500e3 || h < 0
%         pause
%     end
%     
%     wind = atmoshwm(lat,lon,h,'day',site.day,'model','quiet','version','14') * site.wind_mag;
%     uw = wind(1);
%     vw = wind(2);
%     
% else
%     
%     uw = 0;
%     vw = 0;
%     
% end
% 
