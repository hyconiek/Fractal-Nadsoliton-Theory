"""Binary64 interval arithmetic with explicit outward nextafter rounding.

Public endpoints are exact Fractions, so midpoint/Taylor-radius calculations
in the audited source remain exact. Nonfinite arithmetic is rejected.
"""
from fractions import Fraction as F
import math

class Interval:
    def __init__(self,lo,hi=None):
        a=F(lo);b=F(lo if hi is None else hi)
        if a>b:raise ValueError('reversed bounds')
        try:x=float(a);y=float(b)
        except OverflowError as exc:raise ValueError('nonfinite input') from exc
        if not math.isfinite(x) or not math.isfinite(y):raise ValueError('nonfinite input')
        self.a=math.nextafter(x,-math.inf) if F(x)>a else x
        self.b=math.nextafter(y,math.inf) if F(y)<b else y
    @property
    def lo(self):return F(self.a)
    @property
    def hi(self):return F(self.b)
    @staticmethod
    def cast(x):return x if isinstance(x,Interval) else Interval(x)
    @classmethod
    def rounded(cls,a,b):
        if not math.isfinite(a) or not math.isfinite(b):raise ValueError('nonfinite operation')
        o=object.__new__(cls);o.a=math.nextafter(a,-math.inf);o.b=math.nextafter(b,math.inf)
        if not math.isfinite(o.a) or not math.isfinite(o.b):raise ValueError('nonfinite enclosure')
        return o
    def __add__(self,o):
        o=self.cast(o);return self.rounded(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):
        o=object.__new__(Interval);o.a=-self.b;o.b=-self.a;return o
    def __sub__(self,o):return self+-self.cast(o)
    def __rsub__(self,o):return self.cast(o)+-self
    def __mul__(self,o):
        o=self.cast(o);v=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b]
        return self.rounded(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=self.cast(o)
        if o.a<=0<=o.b:raise ValueError('zero denominator')
        v=[self.a/o.a,self.a/o.b,self.b/o.a,self.b/o.b]
        return self.rounded(min(v),max(v))
    def __rtruediv__(self,o):return self.cast(o)/self
    def __pow__(self,n):
        if n<0:return Interval(1)/(self**(-n))
        out=Interval(1)
        for _ in range(n):out=out*self
        return out
