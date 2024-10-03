import numpy as np
from shapely.geometry import Point, LineString, Polygon
from numpy.linalg import norm
from math import pi



def roundness(s):
    '''
    roundness or circularity
    :param s: s should be the list of
    :return:
    '''
    geom = Polygon(s)
    area = geom.area
    c = geom.length
    return 4 * pi * area / np.square(c)

def AngleCal(l1, l2):
    v1 = np.array([l1.xy[0][0] - l1.xy[0][-1], l1.xy[1][0] - l1.xy[1][-1]])
    v2 = np.array([l2.xy[0][0] - l2.xy[0][-1], l2.xy[1][0] - l2.xy[1][-1]])
    angle = (v1[0] * v2[0] + v1[1] * v2[1]) / (norm(v1, ord=2, axis=0) * norm(v2, ord=2, axis=0))
    return angle

def AngleCal_turning(l1, l2, intersctpt):
    '''
    geom of intersection with geoms of relevant links
    :param l1:
    :param l2:
    :param intersctpt:
    :return:
    '''
    t1 = [Point(l1.xy[0][0], l1.xy[1][0]),
          Point(l1.xy[0][-1], l1.xy[1][-1])]
    t2 = [Point(l2.xy[0][0], l2.xy[1][0]),
          Point(l2.xy[0][-1], l2.xy[1][-1])]
    if t1[0].distance(intersctpt) > 1:
        t1.reverse()
    if t2[0].distance(intersctpt) > 1:
        t2.reverse()
    angle = AngleCal(LineString(t1), LineString(t2))
    return angle


def AngleCal_nd(nd1, nd2, nd3):
    v1 = np.array([nd2.x - nd1.x, nd2.y - nd1.y])
    v2 = np.array([nd3.x - nd2.x, nd3.y - nd2.y])
    length = norm(v1, ord=2, axis=0) * norm(v2, ord=2, axis=0)
    if length != 0:
        return (v1[0] * v2[0] + v1[1] * v2[1]) / length
    else:
        raise ValueError("AngleCal_nd problems")