function z_from_time,value
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
redshift = spline(t,z,value)
return,redshift
end
