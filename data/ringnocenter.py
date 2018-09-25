make('Mercury'     ,1.43E+12 ,183       ,1    ,10,10,    '--' )
for i in range(90):
    make('d%s'%i, .8e8 , i*4, 1, 8e6, 1e24, 'Mercury') 
#make('Spaceship', 2.1e8, 0,1,100,100,'Mercury',cls=SpaceShip)
