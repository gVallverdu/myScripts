#!/usr/bin/env python
# -*- coding=utf-8 -*-

""" 
Conversion en coordonnees spheriques
------------------------------------

Relation entre les coordonnées :

.. 

    x = r * sin theta * cos phi
    y = r * sin theta * sin phi
    z = r * cos theta

"""

import doctest
from math import sin, cos, degrees, radians, sqrt, acos, pi

def spheric2cart(r, theta, phi, unit = "deg"):
    """ Calcule les coordonnees cartesiennes x, y, z d'un point M en fonction de ses
    coordonnees spheriques r, theta, phi.

    arguments :

        * r     : rayon
        * theta : angle entre le vecteur OM et l'axe z
        * phi   : angle antre l'axe x et la projection du vecteur OM dans le plan xOy
        * unit  : precise si les angles sont donnes en degres (deg) ou radians (rad)

    """ 

    if unit.strip().lower() == "deg":
        th = radians(theta)
        ph = radians(phi)
    elif unit.strip().lower() == "rad":
        th = theta
        ph = phi
    else:
        raise NameError("unit must take the value 'deg' or 'rad'. You enter " + \
            unit.strip().lower())

    x = r * sin(th) * cos(ph)
    y = r * sin(th) * sin(ph)
    z = r * cos(th)

    return [x, y, z]

def cart2spheric(xyz, unit = "deg", origine = [0., 0., 0.]):
    """ Calcule les coordonnees spheriques r, theta, phi d'un point M en fonction de ses
    coordonnees cartesiennes x, y, z.

    arguments :
        
        * xyz     : coordonnees cartesiennes
        * unit    : precise si les angles sont renvoye en degres (deg) ou radians (rad)
        * origine : coordonnees x, y, z de l'origine du repere

    >>> x =  2. * cos(radians(120.))
    >>> y = -2. * sin(radians(120.))
    >>> z = 0.
    >>> cart2spheric([x,y,z])
    (1.9999999999999998, 90.0, 240.00000000000003)

    """

    # substract origin
    x = xyz[0] - origine[0]
    y = xyz[1] - origine[1]
    z = xyz[2] - origine[2]
     
    # compute coordinates
    r = sqrt(x**2 + y**2 + z**2)
    th = acos(z / r)
    if y >= 0:
        ph = acos(x / r)
    else:
        ph = 2. * pi - acos(x / r)

    # degrees or radians
    if unit.strip().lower() == "rad":
        theta = th
        phi = ph
    elif unit.strip().lower() == "deg":
        theta = degrees(th)
        phi = degrees(ph)
    else:
        raise NameError("unit must take the value 'deg' or 'rad'. You enter " + \
            unit.strip().lower())

    return r, theta, phi

if __name__ == "__main__":
    doctest.testmod()

