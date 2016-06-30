PRO rotatetipsy,filename,angle
;Angle in degrees
angle_rad = angle*!PI/180.0

rtipsy,filename,h,g,d,s
ax = [sin(angle_rad + !PI/2),0,cos(angle_rad + !PI/2)]
ay = [0,1,0]
az = [sin(angle_rad)        ,0,cos(angle_rad)]
basis = [[ax],[ay],[az]]

gpos = [[g.x],[g.y],[g.z]]
dpos = [[d.x],[d.y],[d.z]]
spos = [[s.x],[s.y],[s.z]]

gpos = transpose(transpose(basis)#transpose(gpos))
dpos = transpose(transpose(basis)#transpose(dpos))
spos = transpose(transpose(basis)#transpose(spos))

g.x = gpos[*,0]
g.y = gpos[*,1]
g.z = gpos[*,2]
d.x = dpos[*,0]
d.y = dpos[*,1]
d.z = dpos[*,2]
s.x = spos[*,0]
s.y = spos[*,1]
s.z = spos[*,2]

wtipsy,filename + '.' + strtrim(fix(angle),2),h,g,d,s,/standard

END
