FUNCTION binzread,filename

readcol,filename,i,red,mass,vel,temp,rho,z0,z1,z2,z3,bcoord,bsize,x,y,z,mass4,vel4,rho4,temp4,z4
data = {i:0, mass:0.0, vel:0.0, temp:0.0, rho:0.0, z0:0.0, z1:0.0, z2:0.0, z3:0.0, bcoord:0.0, bsize:0.0, x:0.0, y:0.0, z:0.0, mass4:0.0, vel4:0.0, rho4:0.0, temp4:0.0, z4:0.0}
data = replicate(data, n_elements(i))
data.i    = i
data.mass = mass
data.vel  = vel
data.temp = temp
data.rho  = rho
data.z0   = z0
data.z1   = z1
data.z2   = z2
data.z3   = z3
data.bcoord = bcoord
data.bsize  = bsize
data.x    = x
data.y    = y
data.z    = z
data.mass4 = mass4
data.vel4 = vel4
data.rho4 = rho4
data.temp4 = temp4
data.z4   = z4
RETURN,data
END
