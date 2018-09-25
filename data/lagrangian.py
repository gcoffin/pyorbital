earth=make(name='Earth' ,distance= 1.43E+12 ,phaseindeg= 183 ,direction=1    , radius=6.03E+07,mass=5.69E+27, parentname=   '--' )
r = 3.84E+08
moon=make(name='Moon' ,distance=r2 ,phaseindeg=0 ,direction=1 ,radius =1.74E+06,7.35E+22,parentname='Earth' )
m1 =earth.mass
m2 =moon.mass
r2 = m1/(m1+m2)*r
q = earth.mass/moon.mass
e = (q/3.)**(1./3)
#http://fr.wikipedia.org/wiki/Point_de_Lagrange
l1 = r2 + r*(-e+1./3*e**2+1./9*e**3)

    
