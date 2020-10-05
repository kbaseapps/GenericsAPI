import pandas as pd
import numpy as np



class MatrixValidationException(Exception): pass


def mstr(m):
    maxlen = 400
    mstr = str(m)
    if len(mstr) > maxlen:
        mstr = mstr[:maxlen] + '...'
    return mstr

def rangestr(rng, inclusive):
    s = (
        '%s%g,%g%s'
        % (
            '[' if inclusive[0] else '(',
            rng[0],
            rng[1],
            ']' if inclusive[1] else ')',
        )
    )
    return s

def opstr(opname):
    return ' for operation ' + opname if opname is not None else ''


def assert_is_nonnegative(m, opname=None):
    if isinstance(m, pd.DataFrame):
        m = m.values

    if not (m >= 0).all():
        raise MatrixValidationException(
            'Matrix\n`%s` must have nonnegative entries' # TODO prev stack info?
            % mstr(m)
            + opstr(opname)
        )

def assert_in_range(m, rng, inclusive, opname=None):
    if isinstance(m, pd.DataFrame):
        m = m.values

    if (m < rng[0]).any() or (rng[1] < m).any() or \
       not inclusive[0] and (m == rng[0]).any() or \
       not inclusive[1] and (m == rng[1]).any():
        raise MatrixValidationException(
            'Matrix\n`%s` must have entries in range `%s`' 
            % (mstr(m), rangestr(rng, inclusive))
            + opstr(opname)
        )
