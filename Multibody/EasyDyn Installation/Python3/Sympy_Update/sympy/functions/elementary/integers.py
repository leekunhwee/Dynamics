from __future__ import print_function, division

from sympy.core.basic import C
from sympy.core.singleton import S
from sympy.core.function import Function
from sympy.core import Add
from sympy.core.evalf import get_integer_part, PrecisionExhausted

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################


class RoundFunction(Function):
    """The base class for rounding functions."""

    @classmethod
    def eval(cls, arg):
        if arg.is_integer:
            return arg
        if arg.is_imaginary or (S.ImaginaryUnit*arg).is_real:
            i = C.im(arg)
            if not i.has(S.ImaginaryUnit):
                return cls(i)*S.ImaginaryUnit
            return cls(arg, evaluate=False)

        v = cls._eval_number(arg)
        if v is not None:
            return v

        # Integral, numerical, symbolic part
        ipart = npart = spart = S.Zero

        # Extract integral (or complex integral) terms
        terms = Add.make_args(arg)

        for t in terms:
            if t.is_integer or (t.is_imaginary and C.im(t).is_integer):
                ipart += t
            elif t.has(C.Symbol):
                spart += t
            else:
                npart += t

        if not (npart or spart):
            return ipart

        # Evaluate npart numerically if independent of spart
        if npart and (
            not spart or
            npart.is_real and (spart.is_imaginary or (S.ImaginaryUnit*spart).is_real) or
                npart.is_imaginary and spart.is_real):
            try:
                re, im = get_integer_part(
                    npart, cls._dir, {}, return_ints=True)
                ipart += C.Integer(re) + C.Integer(im)*S.ImaginaryUnit
                npart = S.Zero
            except (PrecisionExhausted, NotImplementedError):
                pass

        spart += npart
        if not spart:
            return ipart
        elif spart.is_imaginary or (S.ImaginaryUnit*spart).is_real:
            return ipart + cls(C.im(spart), evaluate=False)*S.ImaginaryUnit
        else:
            return ipart + cls(spart, evaluate=False)

    def _eval_is_finite(self):
        return self.args[0].is_finite

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real


class floor(RoundFunction):
    """
    Floor is a univariate function which returns the largest integer
    value not greater than its argument. However this implementation
    generalizes floor to complex numbers.

    More information can be found in "Concrete mathematics" by Graham,
    pp. 87 or visit http://mathworld.wolfram.com/FloorFunction.html.

        >>> from sympy import floor, E, I, Float, Rational
        >>> floor(17)
        17
        >>> floor(Rational(23, 10))
        2
        >>> floor(2*E)
        5
        >>> floor(-Float(0.567))
        -1
        >>> floor(-I/2)
        -I

    See Also
    ========

    ceiling
    """
    _dir = -1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            if arg.is_Rational:
                return C.Integer(arg.p // arg.q)
            elif arg.is_Float:
                return C.Integer(int(arg.floor()))
            else:
                return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(C.Integer)[0]

    def _eval_nseries(self, x, n, logx):
        r = self.subs(x, 0)
        args = self.args[0]
        args0 = args.subs(x, 0)
        if args0 == r:
            direction = (args - args0).leadterm(x)[0]
            if direction.is_positive:
                return r
            else:
                return r - 1
        else:
            return r


class ceiling(RoundFunction):
    """
    Ceiling is a univariate function which returns the smallest integer
    value not less than its argument. Ceiling function is generalized
    in this implementation to complex numbers.

    More information can be found in "Concrete mathematics" by Graham,
    pp. 87 or visit http://mathworld.wolfram.com/CeilingFunction.html.

        >>> from sympy import ceiling, E, I, Float, Rational
        >>> ceiling(17)
        17
        >>> ceiling(Rational(23, 10))
        3
        >>> ceiling(2*E)
        6
        >>> ceiling(-Float(0.567))
        0
        >>> ceiling(I/2)
        I

    See Also
    ========

    floor
    """
    _dir = 1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            if arg.is_Rational:
                return -C.Integer(-arg.p // arg.q)
            elif arg.is_Float:
                return C.Integer(int(arg.ceiling()))
            else:
                return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(C.Integer)[1]

    def _eval_nseries(self, x, n, logx):
        r = self.subs(x, 0)
        args = self.args[0]
        args0 = args.subs(x, 0)
        if args0 == r:
            direction = (args - args0).leadterm(x)[0]
            if direction.is_positive:
                return r + 1
            else:
                return r
        else:
            return r
