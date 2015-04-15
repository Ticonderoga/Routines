


def angle_to_point(point, centre):
    '''calculate angle in 2-D between points and x axis'''
    from numpy import arctan,pi

    delta = point - centre
    res = arctan(delta[1] / delta[0])
    if delta[0] < 0:
        res += pi
    return res


def area_of_triangle(p1, p2, p3):
    '''calculate area of any triangle given co-ordinates of the corners'''
    from numpy.linalg import norm
    from numpy import cross

    return norm(cross((p2 - p1), (p3 - p1)))/2.


def convex_hull(points):
    '''Calculate subset of points that make a convex hull around points
    Recursively eliminates points that lie inside two neighbouring points until only convex hull is remaining.
    :Parameters:
        points : ndarray (2 x m)
            array of points for which to find hull
    :Returns:
        hull_points : ndarray (2 x n)
            convex hull surrounding points
    '''
    from numpy import apply_along_axis,asarray
    
    n_pts = points.shape[1]
    centre = points.mean(1)
    angles = apply_along_axis(angle_to_point, 0, points, centre)
    pts_ord = points[:,angles.argsort()]
    pts = [x[0] for x in zip(pts_ord.transpose())]
    prev_pts = len(pts) + 1
    k = 0
    while prev_pts > n_pts:
        prev_pts = n_pts
        n_pts = len(pts)
        i = -2
        while i < (n_pts - 2):
            Aij = area_of_triangle(centre, pts[i],     pts[(i + 1) % n_pts])
            Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], \
                                   pts[(i + 2) % n_pts])
            Aik = area_of_triangle(centre, pts[i],     pts[(i + 2) % n_pts])
            if Aij + Ajk < Aik:
                del pts[i+1]
            i += 1
            n_pts = len(pts)
        k += 1
    return asarray(pts)

#~ if __name__ == "__main__":
    #~ points = random.random_sample((2,40))
    #~ hull_pts = convex_hull(points)
