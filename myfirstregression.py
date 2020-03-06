#myfirstregression.py
import numpy as np
from casadi import *
from scipy.special import hermite
from numpy.polynomial.hermite_e import *
import functools
import operator
import warnings


lis = np.array([ [7,2] , [1,2], [1,2], [7,2], [1,2] ])
  
# using reduce to compute sum of list 
print(functools.reduce(operator.mul,lis))


class RankWarning(UserWarning):
    """Issued by chebfit when the design matrix is rank deficient."""
    pass

print(hermeval(4,2))
print('wut')
np.random.seed(3)
x = DM(np.random.randn(50))
z = DM(np.random.randn(50))
y = DM(np.random.randn(50)+0.001*x)
a = SX.sym('a')
b = SX.sym('b')
c = SX.sym('c', 1, 1)
d = SX.sym('d')
degree = 2
test = SX.sym('test',(degree+1)**2,1)
objective = sum1((y-(a*(x**3-3*x)+b*(x**2-1)+c*x+d))**2)
# objective2 = sum1((y-(hermevander2d(x,z,[degree,degree])@test))**2)
x_0 = DM.ones(4)
# x_0 = DM.ones((degree+1)**2)

print(np.polynomial.hermite_e.HermiteE((1,2,3,4,5)))

nlp = {'x':vertcat(a,b,c,d), 'f':objective}
# solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':0})
# nlp = {'x':test, 'f':objective2}
solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':0})
sol = solver(x0=x_0)['x']
print(sol.shape)
print(sol)
# predictedys = hermevander2d(x,z,[degree,degree])@sol
# print(sum1(fabs(y-predictedys))/5000)
# print(sum(y))
# print(hermevander2d(np.ones(100),np.ones(100),np.array([5,5])).shape)
# print(sol)
# asdf = hermefit([x,z],y,3)
# print(asdf)
jkl = hermevander(DM([1,2,3,4,5,6,7,8]),3)
# print(hermevander2d(np.array([1,2]),np.array([4,5]),np.array([2,2])))
# print(asdf.shape)
# print(jkl.shape)
# print(asdf*jkl)
# np.savetxt("myfirstregression.csv", np.array([x,y,z]), delimiter=",")

print(hermevander(np.ones(2),3))
test = DM([1,2,3,4,5,6,7,8])
print(test[np.array([2,4]),0])
x = hermevander(np.zeros(2),3)
y = x>-0.5
print(x)
jkl = 4
asdf = np.zeros((jkl,jkl),dtype=bool)
for i in range(jkl):
    asdf[i,0:(jkl-i)]=True
print(asdf)
print(asdf.reshape(1,jkl**2))
print(np.arange(jkl**2).reshape(jkl,jkl)[asdf])
print(SX(range(16))[np.arange(jkl**2).reshape(jkl,jkl)[asdf],0])
print(y)
print(x[y])
print(hermevander2d(np.ones(2),np.zeros(2),[3,3]))
print(1 + 'x')

def hermevander_casadiSYM(x, deg):
    
    ideg = operator.index(deg)
    dims = (ideg + 1,) + x.shape
    v = SX.zeros(dims[:2])
    v[0,0:dims[1]] = 1
    if ideg > 0:
        v[1,0:dims[1]] = x.T
        for i in range(2, ideg + 1):
            v[i,0:dims[1]] = (v[i-1,0:dims[1]]*x.T - v[i-2,0:dims[1]]*(i - 1))
    return v.T

asdf = SX.sym('asdf',10,1)
print(hermevander_casadiSYM(asdf,6))
print('x' + 5)

def _nth_slice(i, ndim):
    sl = [np.newaxis] * ndim
    sl[i] = slice(None)
    return tuple(sl)


def _vander_nd(vander_fs, points, degrees):
    
    n_dims = len(vander_fs)
    if n_dims != len(points):
        raise ValueError(
            f"Expected {n_dims} dimensions of sample points, got {len(points)}")
    if n_dims != len(degrees):
        raise ValueError(
            f"Expected {n_dims} dimensions of degrees, got {len(degrees)}")
    if n_dims == 0:
        raise ValueError("Unable to guess a dtype or shape when no points are given")

    # convert to the same shape and type
    # points = tuple(np.array(tuple(points), copy=False) + 0.0)

    # produce the vandermonde matrix for each dimension, placing the last
    # axis of each in an independent trailing axis of the output
    vander_arrays = (
        vander_fs[i](points[i], degrees[i])[(...,) + _nth_slice(i, n_dims)]
        for i in range(n_dims)
    )

    # we checked this wasn't empty already, so no `initial` needed
    return functools.reduce(operator.mul, vander_arrays)


def _vander_nd_flat(vander_fs, points, degrees):
    """
    Like `_vander_nd`, but flattens the last ``len(degrees)`` axes into a single axis
    Used to implement the public ``<type>vander<n>d`` functions.
    """
    v = _vander_nd(vander_fs, points, degrees)
    v = v.reshape(v.shape[:-len(degrees)] + (-1,))
    if len(v.shape)>2.5:
        v = v.reshape(v.shape[0],v.shape[2])
    return v

print('**********')
print(_vander_nd_flat((hermevander,hermevander),[x,z],[5,5]))
# print(hermevander2d(x,z,[5,5]))






def herme2d_fit(x, y, deg, rcond=None, full=False, w=None):
    """
    Helper function used to implement the ``<type>fit`` functions.
    Parameters
    ----------
    vander_f : function(array_like, int) -> ndarray
        The 1d vander function, such as ``polyvander``
    c1, c2 :
        See the ``<type>fit`` functions for more detail
    """
    # x = np.asarray(x) + 0.0
    y = np.asarray(y) + 0.0
    deg = np.asarray(deg)

    # check arguments.
    if deg.ndim > 1 or deg.dtype.kind not in 'iu' or deg.size == 0:
        raise TypeError("deg must be an int or non-empty 1-D array of int")
    if deg.min() < 0:
        raise ValueError("expected deg >= 0")
    # if x.ndim != 1:
        # raise TypeError("expected 1D vector for x")
    # if x.size == 0:
    #     raise TypeError("expected non-empty vector for x")
    if y.ndim < 1 or y.ndim > 2:
        raise TypeError("expected 1D or 2D array for y")
    # if len(x) != len(y):
        # raise TypeError("expected x and y to have same length")

    # if deg.ndim == 0:
    #     lmax = deg
    #     order = lmax + 1
    #     van = vander_f(x, lmax)

    # else:
    #     deg = np.sort(deg)
    #     lmax = deg[-1]
    #     order = len(deg)
    #     van = vander_f(x, lmax)[:, deg]
    van = _vander_nd_flat((hermevander,hermevander),x,[deg,deg])
    order = (deg+1)**2


    # set up the least squares matrices in transposed form
    lhs = van.T
    rhs = y.T
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            raise TypeError("expected 1D vector for w")
        if len(x) != len(w):
            raise TypeError("expected x and w to have same length")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    if rcond is None:
        rcond = len(y)*np.finfo(y.dtype).eps

    # Determine the norms of the design matrix columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1
    lhsT = lhs.T
    if len(lhsT.shape)>2.5:
        lhsT = lhsT.reshape(lhsT.shape[0],lhsT.shape[2])
        scl = np.sqrt(np.sum(np.square(lhs),axis=2)).T
        lhsT_over_scl = lhsT/scl
    else:
        lhsT_over_scl = lhsT/scl

    # Solve the least squares problem.
    c, resids, rank, s = np.linalg.lstsq(lhsT_over_scl, rhs.T, rcond)
    c = (c.T/scl).T

    # Expand c to include non-fitted coefficients which are set to zero
    if deg.ndim > 0:
        if c.ndim == 2:
            cc = np.zeros((lmax+1, c.shape[1]), dtype=c.dtype)
        else:
            cc = np.zeros(lmax+1, dtype=c.dtype)
        cc[deg] = c
        c = cc

    # warn on rank reduction
    if rank != order and not full:
        msg = "The fit may be poorly conditioned"
        warnings.warn(msg, RankWarning, stacklevel=2)

    if full:
        return c, [resids, rank, s, rcond]
    else:
        return c



degree = 2
coefs = herme2d_fit([x,z],y,degree).T
print(coefs.shape)
if coefs.shape[0] > 1:
    coefs = coefs.T
vals = _vander_nd_flat((hermevander,hermevander),[x,z],[degree,degree])
print((coefs*vals).shape)

