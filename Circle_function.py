from numpy import *

def distance(x1, y1, x2, y2):
    '''
    This function takes in two coordinates and returns the Euclidean
    distance between them

    Parameters
    ----------

    x1 : int
        first x coordinate
    y1 : int
        first y coordinate
     x2 : int
        second x coordinate
    y2 : int
        second y coordinate
    
    Returns
    -------
        This function returns the distance between the two coordinates.
    '''

    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

def circle_spiral(params, length, percentage, segment, Vm, phi, theta):

    '''
    This function maps impurities to the perimeter of a circle/arc.
    The spin orientation is helical.
    
    Parameters
    ----------
    params : dict
        dictionary containing system profile
    length : int
        radius or arc length
    percentage : int
        percentage of circle removed (0% -> circle, 50% -> semi-circle)
    segment : str
        'radius' or 'arc'
    Vm : float
        magnetic scattering strength of impurities
    phi : float
        azimuthal angle
    theta : float
        polar angle
        
    Returns
    -------
        array of characteristics for each impurity
        
    Some of this code was first developed by Dr. Bernd Braunecker but has been
    heavily modified for use in this project
    '''    
    gap = pi * (percentage/100)

    x = -1
    y = -1

    Nx2 = params['Nx']//2
    Ny2 = params['Ny']//2

    impurities = []

    i = 1

    # optimal pitch shift
    #dphi = 2*arccos((-2*params['t'] - params['mu'])/(2*params['t']))
    dphi=2*params["km"]

    if gap ==0:
        # approximation of optimal pitch shift
        INT = round(dphi * 320/(2*pi))
        dphi = (2*pi*INT)/320

    if segment == 'radius':
        # create impurities in a circular structure
        r  = length
    elif segment == 'arc':
        # create impurities in a arc structure
        if gap >= 0 and gap < pi:
            r = length/((pi - gap)*2)
        else:
            for idx in range(length):
                        yn = (Ny2 - length//2) + idx
                        impurities.append([Nx2, yn, Vm, idx*dphi, theta, 0.])
            return impurities

    # circle becomes symmetric if number of subdivisions of ph
    # is of the form (4*n)

    for ph in linspace(0+gap,2*pi-gap,400,endpoint=True):

        # determines the offset that centres the chain of impurities
        # for gap of pi, no impurities are placed

        length = len(impurities)

        if (gap >= 0) and  (gap < pi/2):
            xn = Nx2 + round((1/2)*r*(1-cos(gap))) + round(r*cos(ph))
            yn = Ny2 + round(r*sin(ph))
        elif (gap >= pi/2) and (gap < pi):
            xn = Nx2 + round((1/2)*r*(1 + sin(gap - pi/2))) + round(r*cos(ph))
            yn = Ny2 + round(r*sin(ph))
        else:
            break

        if ph == 0+gap:
            x = xn
            y = yn
            impurities.append([x, y, Vm, phi, theta, 0.])
            continue

        # test if new position has already been placed and if not, determines
        # the distance between the new impurity and the previously placed one
        # to determine if there is a diagonal gap.
        # if not, an impurity is placed. 
        if ((xn!=x) or (yn!=y)) and (distance(xn, yn, x, y) == 1):
            x = xn
            y = yn
            impurities.append([x,y,Vm,i*dphi,theta,0.])
            i+=1
        # if there is a diagonal, we first determine where in the circle the
        # the impurity is being placed which then tells us which direction to
        # place the new impurity such that it is on the inside of the circle.
        elif ((xn!=x) or (yn!=y)) and (distance(xn, yn, x, y) > 1):
            if (ph >= 0) and (ph < pi/2):
                impurities.append([x - 1, y, Vm, i*dphi, theta, 0.])
            elif (ph >= pi/2) and (ph < pi):
                impurities.append([x, y - 1, Vm, i*dphi, theta, 0.])
            elif (ph >= pi) and (ph < 3*pi/2):
                impurities.append([x + 1, y, Vm, i*dphi, theta, 0.])
            elif (ph >= 3*pi/2) and (ph < 2*pi):
                impurities.append([x, y + 1, Vm, i*dphi, theta, 0.])
            i+=1

            x = xn
            y = yn
            impurities.append([x,y,Vm,i*dphi,theta,0.])
            i+=1
    if gap == 0:
        impurities.pop()

    return impurities