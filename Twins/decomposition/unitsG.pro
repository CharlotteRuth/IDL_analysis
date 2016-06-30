pro unitsG,gas,dark,stars,a

ms_units=13.598e16                 ;2304Gal1
dist_units=1e5*a                  
v_units=2417.0*a


stars.mass=ms_units*stars.mass
gas.mass=ms_units*gas.mass
dark.mass=ms_units*dark.mass

stars.tform=13.7/0.33*stars.tform
stars.phi=v_units*v_units*stars.phi

gas.x=dist_units*gas.x
gas.y=dist_units*gas.y
gas.z=dist_units*gas.z
dark.x=dist_units*dark.x
dark.y=dist_units*dark.y
dark.z=dist_units*dark.z
stars.x=dist_units*stars.x
stars.y=dist_units*stars.y
stars.z=dist_units*stars.z

gas.vx=v_units*gas.vx
gas.vy=v_units*gas.vy
gas.vz=v_units*gas.vz
dark.vx=v_units*dark.vx
dark.vy=v_units*dark.vy
dark.vz=v_units*dark.vz
stars.vx=v_units*stars.vx
stars.vy=v_units*stars.vy
stars.vz=v_units*stars.vz

return

end
