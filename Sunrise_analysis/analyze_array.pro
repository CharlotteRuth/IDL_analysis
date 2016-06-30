particles = [-1]
for i = 0L, 2219628 do $
  if (where(data[i].L_LAMBDA lt 0))[0] ne -1 then particles = [particles ,i] 
