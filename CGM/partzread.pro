FUNCTION partzread,filename

readcol,filename,mass,rho,temp,z0,z1,z2,hsmooth,x,y,z,vx,vy,vz
data = {mass:0.0, rho:0.0, temp:0.0, z0:0.0, z1:0.0, z2:0.0, hsmooth:0.0, x:0.0, y:0.0, z:0.0, vx:0.0, vy:0.0, vz:0.0}
data = replicate(data,n_elements(mass))
data.mass = mass
data.rho  = rho
data.temp = temp
data.z0   = z0
data.z1   = z1
data.z2   = z2
data.hsmooth = hsmooth
data.x    = x
data.y    = y
data.z    = z
data.vx   = vx
data.vy   = vy
data.vz   = vz

RETURN,data
END
