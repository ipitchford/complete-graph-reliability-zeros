"""Run the published phase tests in a second arithmetic: dyadic integer intervals.
Only integer arithmetic is used after parsing decimal inputs as exact rationals.
This checks arithmetic, retaining the certificate's geometry and subdivision.
"""
from fractions import Fraction
from functools import lru_cache
from math import isqrt, factorial
from pathlib import Path
import sys


class CertificateFailure(ArithmeticError):
    """Raised when a certified inequality fails; never disabled by python -O."""


def require(condition, message='certificate inequality failed'):
    if not condition:
        raise CertificateFailure(message)


Q=1<<220
def ceildiv(a,b): return -((-a)//b)
class R:
    def __init__(self,v=0):
        if isinstance(v,R): self.lo,self.hi=v.lo,v.hi
        elif isinstance(v,(list,tuple)): self.lo,self.hi=R(v[0]).lo,R(v[1]).hi
        else:
            v=Fraction(v)
            self.lo=v.numerator*Q//v.denominator
            self.hi=ceildiv(v.numerator*Q,v.denominator)
    @classmethod
    def raw(cls,lo,hi):
        require(lo<=hi)
        v=object.__new__(cls);v.lo,v.hi=lo,hi
        return v
    def __add__(a,b):
        if isinstance(b,(Z,complex)):return NotImplemented
        b=R(b);return R.raw(a.lo+b.lo,a.hi+b.hi)
    __radd__=__add__
    def __neg__(a):return R.raw(-a.hi,-a.lo)
    def __sub__(a,b):return a+(-R(b))
    def __rsub__(a,b):return R(b)+(-a)
    def __mul__(a,b):
        if isinstance(b,(Z,complex)):return NotImplemented
        b=R(b);t=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi]
        return R.raw(min(t)//Q,ceildiv(max(t),Q))
    __rmul__=__mul__
    def __truediv__(a,b):
        b=R(b);require(b.lo*b.hi>0)
        return a*R.raw(Q*Q//b.hi,ceildiv(Q*Q,b.lo))
    def __rtruediv__(a,b):return R(b)/a
    def __pow__(a,k):
        require(k>=0)
        if k==0:return R(1)
        if k%2==0:
            t=(a**(k//2))
            x,y=t.lo,t.hi
            return R.raw((0 if x*y<=0 else min(x*x,y*y))//Q, ceildiv(max(x*x,y*y),Q))
        return a*(a**(k-1))
    def __abs__(a):return R.raw(0 if a.lo*a.hi<=0 else min(abs(a.lo),abs(a.hi)),max(abs(a.lo),abs(a.hi)))
    def __gt__(a,b):return a.lo>R(b).hi
    def __lt__(a,b):return a.hi<R(b).lo
    def __contains__(a,b):return a.lo<=b.lo and b.hi<=a.hi
    @property
    def a(self):return self.lo/Q
    def __repr__(a):return f'[{a.lo/Q:.8g}, {a.hi/Q:.8g}]'
def sqrt(a):
    a=R(a); require(a.lo>=0)
    return R.raw(isqrt(a.lo*Q),isqrt(a.hi*Q)+1)
class Z:
    def __init__(self,v=0,y=0):
        if isinstance(v,Z):self.real,self.imag=v.real,v.imag
        elif isinstance(v,complex): self.real,self.imag=R(v.real),R(v.imag)
        else:self.real,self.imag=R(v),R(y)
    def __add__(a,b):
        b=Z(b);return Z(a.real+b.real,a.imag+b.imag)
    __radd__=__add__
    def __neg__(a):return Z(-a.real,-a.imag)
    def __sub__(a,b):return a+-Z(b)
    def __rsub__(a,b):return Z(b)+-a
    def __mul__(a,b):
        b=Z(b);return Z(a.real*b.real-a.imag*b.imag,a.real*b.imag+a.imag*b.real)
    __rmul__=__mul__
    def __truediv__(a,b):
        b=Z(b);u=a*Z(b.real,-b.imag);d=b.real**2+b.imag**2
        return Z(u.real/d,u.imag/d)
    def __pow__(a,k):
        if k==0:return Z(1)
        if k==1:return a
        if k%2:return a*a**(k-1)
        u=a**(k//2); return Z(u.real**2-u.imag**2,2*u.real*u.imag)
    def __abs__(a):return sqrt(a.real**2+a.imag**2)
def atan_recip(k):
    total=R(0)
    for j in range(100):
        total+=R(Fraction((-1)**j,(2*j+1)*k**(2*j+1)))
    err=R(Fraction(1,201*k**201))
    return total+R.raw(-err.hi,err.hi)
pi=16*atan_recip(5)-4*atan_recip(239)
hp=pi/2
@lru_cache(maxsize=100000)
def sincos_point(x):
    k=(2*x+(hp.lo+hp.hi)//2)//(hp.lo+hp.hi)
    u=R.raw(x,x)-k*hp
    require(abs(u)<1)
    t,c=R(1),R(1)
    v,s=u,u
    step=-(u**2)
    for j in range(1,28):
        t=t*step/((2*j-1)*2*j); c+=t
        v=v*step/(2*j*(2*j+1)); s+=v
    err=R(Fraction(1,factorial(55)))
    both=R.raw(-err.hi,err.hi);c+=both;s+=both
    return [(s,c),(c,-s),(-s,-c),(-c,s)][k%4]
def trigo(a,index):
    a=R(a);ends=[sincos_point(x)[index] for x in (a.lo,a.hi)];lo=min(t.lo for t in ends);hi=max(t.hi for t in ends)
    ks=[2*x//(hp.lo+hp.hi) for x in (a.lo,a.hi)]
    for k in range(ks[0]-1,ks[1]+3):
        p=k*hp
        if p.hi<a.lo or a.hi<p.lo:continue
        extrem=(k+index)%4
        if extrem==1:hi=Q
        if extrem==3:lo=-Q
    return R.raw(lo,hi)
class Arithmetic:
    dps=0
    mpf=R
    mpc=Z
    pi=pi
    sqrt=staticmethod(sqrt)
    sin=staticmethod(lambda x:trigo(x,0))
    cos=staticmethod(lambda x:trigo(x,1))
    @staticmethod
    def exp(x):
        require(x==1)
        return R([2,3])
iv=Arithmetic()

if __name__ == "__main__":
    import runpy
    from types import ModuleType

    if sys.argv[1:] not in ([], ["--arc"]):
        raise SystemExit("Usage: python3 check_intervals_independent.py [--arc]")
    # Supply the separate arithmetic implementation to the same geometric tests.
    backend = ModuleType("certified_intervals")
    backend.iv = iv
    sys.modules["certified_intervals"] = backend
    script = "check_arc.py" if sys.argv[1:] else "whittaker_certificate.py"
    runpy.run_path(str(Path(__file__).with_name(script)), run_name="__main__")
