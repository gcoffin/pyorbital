This program is a simulator of gravity for a system of N bodies as well as a space trip simulator.

How to install
--------------

You will need 

- python: this software has only been tested with python 2.7. I
  haven't tested it with python 3 but apparently that shouldn't be too
  much work to adapt.

- numpy: https://pypi.python.org/pypi/numpy (a reasonably recent version should work)

- pygame: http://www.pygame.org/download.shtml

How to run
----------
python orbital.py
or
python orbital.py data\<datafile>

The datafile is a csv or python file that contains the configuration
(position, size and mass of the bodies, including the spaceship).
There are examples containing the solar system as well as some random
generation.

How to play
----------
- up arrow to thrust

- down arrow to brake

- left and right to turn

- Because flight in outer space can be very long, you can increase the
  simulation speed. Increase with l (L key) and decrease with k.  The
  hud will show the time acceleration factor. Above a certain level
  the system will have trouble to run the computation and the
  animation will appear less fluid.

- the thrust power can be increased with b and decreased with
  v. You'll need large thrust when taking off but normally low thrust
  will do. Beware also to adapt the time acceleration.

- there is an autozoom feature that is toggled with the key a. If
  centers the screen on either the rocket, the closest planet or a
  point on the surface of that planet. It also adapts the scale to see
  the zone of interest and tries to adapt the time acceleration factor
  to keep things interesting. If disabled, will follow the spacecraft.

Add more images
---------------

The system will look for a file in the directory images with a name
such as <name of the planet><diameter>.<extension>, where

- name of the planet is the same name as in the configuration file

- diameter is the apparent diameter of the planet in pixel. If the
image contains the planet such as the circle goes from left to righ
and from bottom to top, you can omit that number. Otherwise (for
instance if there's an atmosphere or rings), use the apparent diameter
on the image.  The planet must be centered.
