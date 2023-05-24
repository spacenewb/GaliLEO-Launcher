function [uw,vw] = wind_generator_para(z,Y,site)

if z < 0
    z = 0;
    
end

if z < 50000
% diifferentiate for para simulation 
[lat,lon,h] = ned2geodetic(Y(1), Y(2), -z, site.lat0, site.lon0, -site.z0, wgs84Ellipsoid);
wind = atmoshwm(lat,lon,h,'day',site.day,'model','quiet','version','14') * site.wind_mag;
uw = wind(1);
vw = wind(2);

else 
    
    uw = 0;
    vw = 0;
    
end
