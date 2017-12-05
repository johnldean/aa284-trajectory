import poliastro as pa
import poliastro.twobody
import poliastro.plotting
import poliastro.maneuver
from poliastro.bodies import Mars 
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt

START_R = 210 
START_INC = 0
M_init = 56
M = M_init
T = 0.014
Isp = 1375
g0 = 9.81


r = [1,0,0]* (START_R*u.km+Mars.R)
v1 = (Mars.k/(Mars.R + START_R*u.km))**0.5
print(v1)
v = [0,1,0]* v1 
v = v.to(u.km/u.s)
ss = pa.twobody.Orbit.circular(Mars,alt=START_R*u.km,inc=30*u.deg)

timestep = 0.01
t = 0
Dv_total = 0
t_ = []
x_ = np.array(ss.r)[np.newaxis]

while (ss.a - Mars.R) < 600 * u.km:
	ve = ss.v/np.linalg.norm(ss.v)/1000
	Dv = ve * T/M*(timestep*60*60*24)
	man = poliastro.maneuver.Maneuver.impulse(Dv)
	ss = ss.apply_maneuver(man)
	ss = ss.propagate(timestep*u.day)
	print(ss.a - Mars.R,M)

	Dv_total += T/M*(timestep*60*60*24)
	M -= T*(timestep*60*60*24)/(Isp*9.81)
	t += timestep
	x_ = np.append(x_,ss.r[np.newaxis],axis=0)
	t_.append(t)


while(ss.inc > 1e-2 * u.rad):
	ve = ss.v/np.linalg.norm(ss.v)
	re = ss.r/np.linalg.norm(ss.r)
	vxr = np.cross(ve,re)	
	if np.sign(vxr[2]) == np.sign(ss.v[2]):
		vxr = -vxr
	vxr = vxr * u.m / u.s
	Dv = vxr * T/M*(timestep*60*60*24)
	man = poliastro.maneuver.Maneuver.impulse(Dv)
	ss = ss.apply_maneuver(man)
	ss = ss.propagate(timestep*u.day)
	print(ss.inc,M)

	Dv_total += T/M*(timestep*60*60*24)
	M -= T*(timestep*60*60*24)/(Isp*9.81)
	t += timestep
	x_ = np.append(x_,ss.r[np.newaxis],axis=0)
	t_.append(t)

print(t,Dv_total,M_init - M)
print(np.log(M_init/M)*Isp * 9.81)
print(g0*Isp*(M_init - M)/T/(24*60*60))
'''
op = pa.plotting.OrbitPlotter()
print(ss.ecc)
#op.plot(ss)
plt.plot(x_[:,0],x_[:,1])
plt.show()
'''