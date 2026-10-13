"""Generates the base sine curve that course 0 of the wall follows.

This is the ONLY place the curve is generated; ``brick_wall_populate`` consumes
it via its ``base_curve`` input.  The brick count defines the wall's arc length:
given the brick length and gap, the target arc length is
``bricks_per_course * (brickL + brick_gap)``, and the sine's Y-extent is scaled
so the curve's ARC length equals it (accounting for the given amplitude, i.e.
a more curved wave needs a shorter Y-extent).

Inputs:
    bricks_per_course: Number of bricks per course (defines the arc length)
    brickL:            Brick length along the wall run [mm]
    brick_gap:         Gap between bricks [mm]
    amplitude:         Sine wave amplitude [mm]
    cycles:            Number of sine cycles across the wall
    samples:           Curve sample count
Outputs:
    Curve:       Rhino sine curve whose arc length = bricks_per_course*(brickL+brick_gap)
    ArcLength:   The resulting arc length [mm]
"""

__author__ = "aa9663"
__version__ = "2023.09.11"

import math
import Rhino.Geometry as rg

bricks    = globals().get("bricks_per_course", 10)
brickL    = globals().get("brickL", 190.0)        # mm
brick_gap = globals().get("brick_gap", 10.0)      # mm
amp       = globals().get("amplitude", 40.0)      # mm
cycles    = globals().get("cycles", 1.0)
samples   = globals().get("samples", 200)

# The brick count defines the wall's ARC length.
target_arc = bricks * (brickL + brick_gap)

def _sine(y_extent):
    pts = []
    for i in range(samples + 1):
        t = i / float(samples)
        x = amp * math.sin(2.0 * math.pi * cycles * t)
        y = t * y_extent
        pts.append(rg.Point3d(x, y, 0.0))
    return rg.Curve.CreateInterpolatedCurve(pts, 3)

# Scale the Y-extent so the curved (arc) length equals target_arc, for this
# amplitude.  The more curved the sine, the shorter the Y-extent needs to be.
lo = max(target_arc * 0.05, 1e-6)
hi = target_arc * 10.0
for _ in range(50):
    mid = 0.5 * (lo + hi)
    if _sine(mid).GetLength() < target_arc:
        lo = mid
    else:
        hi = mid

Curve = _sine(0.5 * (lo + hi))
ArcLength = Curve.GetLength()      # ~= bricks * (brickL + brick_gap)