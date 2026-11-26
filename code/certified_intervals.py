"""Small outward-rounded integer interval backend for the W_0 certificate.

An interval (lo,hi) represents [lo/2^176, hi/2^176]. All rounding is directed.
The only transcendental primitives are pi from Machin's alternating series
and sin/cos from Taylor's theorem. Python integers have unbounded precision.
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial, isqrt


class CertificateFailure(ArithmeticError):
    """Raised when a certified inequality fails; never disabled by python -O."""


def require(condition, message='certificate inequality failed'):
    if not condition:
        raise CertificateFailure(message)


SCALE = 1 << 176


def upper_quotient(a, b):
    return -((-a) // b)


class Real:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi):
        require(lo <= hi)
        self.lo, self.hi = lo, hi

    def __add__(self, other):
        if isinstance(other, (Complex, complex)):
            return NotImplemented
        b = real(other)
        return Real(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Real(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, (Complex, complex)):
            return NotImplemented
        b = real(other)
        corners = [x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return Real(min(corners)//SCALE, upper_quotient(max(corners), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (Complex, complex)):
            return NotImplemented
        b = real(other)
        require(b.lo > 0 or b.hi < 0)
        quotients = [(x*SCALE, y) for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return Real(min(u//v for u,v in quotients), max(upper_quotient(u,v) for u,v in quotients))

    def __rtruediv__(self, other):
        return real(other)/self

    def __pow__(self, k):
        require(isinstance(k,int) and k >= 0)
        if not k:
            return real(1)
        if k % 2:
            lo, hi = self.lo**k, self.hi**k
        else:
            a = abs(self)
            lo, hi = a.lo**k, a.hi**k
        scale = SCALE**(k-1)
        return Real(lo//scale, upper_quotient(hi, scale))

    def __abs__(self):
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return Real(0,max(-self.lo,self.hi))

    def __lt__(self, other):
        return self.hi < real(other).lo

    def __gt__(self, other):
        return self.lo > real(other).hi

    def __contains__(self, other):
        b = real(other)
        return self.lo <= b.lo and b.hi <= self.hi

    @property
    def a(self):
        # Display a certified lower endpoint, rounded down to 15 decimal places.
        scaled = (self.lo * 10**15) // SCALE
        sign = '-' if scaled < 0 else ''
        return f'{sign}{abs(scaled)//10**15}.{abs(scaled)%10**15:015d}'


def real(value):
    if isinstance(value, Real):
        return value
    if isinstance(value, (tuple, list)):
        return Real(real(value[0]).lo, real(value[1]).hi)
    if isinstance(value, int):
        return Real(value*SCALE,value*SCALE)
    x = Fraction(value)*SCALE
    return Real(x.numerator//x.denominator, upper_quotient(x.numerator,x.denominator))


def sqrt(x):
    x = real(x)
    require(x.lo >= 0)
    low, hi = isqrt(x.lo*SCALE), isqrt(x.hi*SCALE)
    return Real(low,hi+(hi*hi<x.hi*SCALE))


class Complex:
    __slots__ = ('real','imag')

    def __init__(self, x=0, y=0):
        if isinstance(x, (complex,Complex)):
            self.real, self.imag = real(x.real),real(x.imag)
        else:
            self.real, self.imag = real(x),real(y)

    def __neg__(self):
        return Complex(-self.real,-self.imag)

    def __add__(self, other):
        b = Complex(other)
        return Complex(self.real+b.real,self.imag+b.imag)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        b = Complex(other)
        return Complex(self.real*b.real-self.imag*b.imag,self.real*b.imag+self.imag*b.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (complex,Complex)):
            b = Complex(other)
            norm_squared = b.real**2+b.imag**2
            return (self*Complex(b.real,-b.imag))/norm_squared
        return Complex(self.real/other,self.imag/other)

    def __pow__(self, k):
        require(isinstance(k,int) and k>=0)
        result, base = Complex(1), self
        while k:
            if k % 2:
                result = result*base
            k //= 2
            if k:
                base = base*base
        return result

    def __abs__(self):
        return sqrt(self.real**2+self.imag**2)


def atan_inverse(m):
    # The even partial sum (through a negative term) is a lower bound.
    terms = 60
    lower = sum(Fraction((-1)**j, (2*j+1)*m**(2*j+1)) for j in range(terms))
    upper = lower + Fraction(1, (2*terms+1)*m**(2*terms+1))
    return real([lower,upper])


PI = 16*atan_inverse(5)-4*atan_inverse(239)


@lru_cache(maxsize=8192)
def trig_point(mid):
    # Quarter turns reduce the argument to [-pi/4,pi/4]. Any integer k gives
    # the exact identity used here; the check independently verifies |u|<1.
    pi_mid = (PI.lo+PI.hi)//2
    k = (4*mid+pi_mid)//(2*pi_mid)
    u = Real(mid,mid) - k*PI/2
    require(abs(u)<1)
    u2 = u**2
    sine_term, cosine_term = u, real(1)
    sine, cosine = sine_term, cosine_term
    for j in range(1,20):
        sine_term = -sine_term*u2/((2*j)*(2*j+1))
        cosine_term = -cosine_term*u2/((2*j-1)*(2*j))
        sine = sine+sine_term
        cosine = cosine+cosine_term
    # Degrees 39 and 38; all real trig derivatives have absolute value <=1.
    err = real([-1,1])/factorial(39)
    sine, cosine = sine+err, cosine+err
    return ((sine,cosine),(cosine,-sine),(-sine,-cosine),(-cosine,sine))[k%4]


def trig(x, component):
    x = real(x)
    mid = (x.lo+x.hi)//2
    h = max(x.hi-mid,mid-x.lo)
    # A first-derivative bound of 1 covers uncertainty in the argument.
    return trig_point(mid)[component] + Real(-h,h)


class Interface:
    mpf = staticmethod(real)
    mpc = Complex
    sqrt = staticmethod(sqrt)
    sin = staticmethod(lambda x: trig(x,0))
    cos = staticmethod(lambda x: trig(x,1))
    pi = PI


iv = Interface()
