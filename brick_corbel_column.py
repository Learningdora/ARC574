"""Optimises a single corbelled column of bricks to maximise lean angle.

The brick's length (``brickL``) runs along the wall (the running/bond direction),
so the wall thickness and the corbel overhang direction use the brick width
(``brickW``).  The odd courses are offset by half the brick length along the wall
run to form a running bond, and that reduced bearing width is fed back into the
stability check.

The column leans outwards: every brick overhangs the brick below it by an amount
``delta`` in the overhang direction, so the column forms a straight line at angle
``theta = atan(delta / (brick_height + mortar_thickness))`` from vertical.

The column is stable while BOTH of these hold at every bed joint k
(the interface carrying the bricks at and above course k):

  1. Strength  - the net tensile stress at the heel of the joint does not
     exceed the mortar's flexural tensile (bond) strength
       M = P * e                       (overturning moment, P = weight above)
       Z = b * d^2/6 ,  A = b * d      (section modulus, joint area)
       sigma_t = M / Z - P / A <= f_t  (the weight's compression offsets the
                                        bending tension; no tension while
                                        e <= d/6)
  2. Stability - the resultant of the weight above the joint stays over the
     supporting brick below (no overturning):
       e <= d / 2
     plus the bricks must still overlap (delta < brick width).

Additionally the whole corbelled wall thickness is capped:
       thickness = brick_width + (n_courses - 1) * delta
                  <= max_thickness_ratio * brick_width   (default 3).

The optimisation maximises delta (== the lean angle) subject to all three
limits; whichever binds first governs the maximum lean.

Module-level functions are pure Python (testable outside Rhino); the
``__main__`` block maps Grasshopper/rhinoscriptsyntax inputs and outputs.

Inputs (brick dimensions and mortar thickness are in mm, converted to metres
internally for the SI solver):
    brickL:          Brick length - runs along the wall (bond direction) [mm]
    brickW:          Brick width - the wall thickness / overhang direction [mm]
    brickH:          Brick height [mm]
    mortar_tensile:  Flexural tensile (bond) strength of the mortar joint [Pa],
                     default 0.2e6 (typical range 0.1-0.4 MPa)
    n_courses:       Number of brick courses in the column (int)
    max_thickness_ratio:  Wall thickness cap in brick widths (default 3)
    brick_density:   Brick density [kg/m3], default 700
    mortar_t:        Mortar bed joint thickness [mm], default 10
    max_overhang:    Cap on each brick's overhang in mm (default = brick width)
    bond_offset:     Lateral offset (mm) of the odd courses along the wall run
                     (Y axis) to make a running bond - half the brick length by
                     default, kept under the brick length so the courses still
                     overlap.  The reduced bearing width is used in the check
    profile:         0 = uniform overhang (straight lean), 1 = curved corbel
                     (ogee profile: small overhang near the base, larger near
                     the top, so the column leans farther) (default 0)
Outputs (geometry and numeric lengths are in mm, matching the input model
unit; angles and stress ratios are dimensionless):
    LeanAngle:       Maximised column lean angle from vertical [deg]
    Overhangs:       Per-course overhang offsets (course 0 = 0.0) [mm]
    PerBrickX:       Horizontal (overhang) centre of each course [mm]
    Frames:          Rhino centre frame (plane) of each brick [mm]; the odd
                     courses are shifted by bond_offset along Y (running bond)
    StressRatio:     Governing net tensile joint stress / flexural tensile
                     strength (<= 1; 0 when every joint is fully compressed)
    GoverningJoint:  Index (course) of the governing joint
    Ok:              True if a stable lean was found
"""

__author__ = "aa9663"
__version__ = "2023.09.11"

import math

G = 9.81
BRICK_DENSITY = 700.0
DEFAULT_MORTAR_T = 0.010
M_TO_MM = 1000.0


def brick_weight(L, W, H, density=BRICK_DENSITY, g=G):
    """Weight (N) of a single brick -> density * volume * g."""
    return density * L * W * H * g


def joint_stress(P, ecc, contact, width):
    """Net tensile stress (Pa) at the heel of a bed joint.

    Bending tension M/Z minus the compression P/A from the weight above,
    clamped at 0 when the whole joint stays in compression (|e| <= d/6).
    """
    if contact <= 0.0 or width <= 0.0:
        return float("inf")
    A = width * contact
    Z = width * contact * contact / 6.0
    return max(P * abs(ecc) / Z - P / A, 0.0)


def course_centres(delta, N):
    """Overhang-direction centre (x) of each course for a uniform overhang."""
    return [i * delta for i in range(N)]


def joint_results(delta, L, W, H, N, density=BRICK_DENSITY, g=G, bond_offset=0.0):
    """Evaluate every bed joint for a uniform overhang delta.

    Joint k carries the bricks at and above course k (k = 0 is the base joint
    under the first course, which has full contact depth W).  Returns a list
    of dicts with P (weight above), e (eccentricity of its resultant from the
    joint centre), contact (overlap depth), moment, Z (section modulus) and the
    net tensile stress.

    ``bond_offset`` shifts the odd courses along the wall face (Y) for a running
    bond, reducing the effective joint width (and hence the section modulus) at
    every intermediate joint by that amount.
    """
    Wb = brick_weight(L, W, H, density, g)
    w_int = max(L - bond_offset, 0.0)     # bearing width at running-bond joints
    res = []
    # Base joint under course 0 (full contact, no overhang on the host).
    contact = W
    ecc = delta * (N - 1) / 2.0
    P = N * Wb
    Z = L * contact * contact / 6.0
    res.append(
        dict(k=0, P=P, ecc=ecc, contact=contact,
             moment=P * ecc, Z=Z,
             stress=joint_stress(P, ecc, contact, L))
    )
    # Intermediate joints under course k (contact depth reduced by overhang).
    for k in range(1, N):
        n_above = N - k
        contact = W - delta
        ecc = delta * n_above / 2.0
        P = n_above * Wb
        Z = w_int * contact * contact / 6.0 if contact > 0 else 0.0
        res.append(
            dict(k=k, P=P, ecc=ecc, contact=contact,
                 moment=P * ecc, Z=Z,
                 stress=joint_stress(P, ecc, contact, w_int))
        )
    return res


def delta_feasible(delta, L, W, H, N, f_t, density, g, max_thickness_ratio,
                   max_overhang=None, bond_offset=0.0):
    """True if a uniform overhang delta keeps the column stable and within cap."""
    if delta < 0.0:
        return False
    if delta >= W:                      # bricks no longer overlap
        return False
    if max_overhang is not None and delta > max_overhang:
        return False
    if N < 1:
        return False
    joints = joint_results(delta, L, W, H, N, density, g, bond_offset=bond_offset)
    # Strength: net joint tension must not exceed flexural tensile strength.
    if max(j["stress"] for j in joints) > f_t:
        return False
    # Stability: resultant of the weight above stays over the supporting brick.
    for j in joints:
        if j["ecc"] > j["contact"] / 2.0:
            return False
    # Wall thickness cap.
    thickness = W + (N - 1) * delta
    if thickness > max_thickness_ratio * W:
        return False
    return True


def max_uniform_overhang(L, W, H, N, f_t, density=BRICK_DENSITY, g=G,
                         max_thickness_ratio=3.0, max_overhang=None,
                         bond_offset=0.0, tol=1e-7):
    """Binary-search the largest uniform overhang that is still stable."""
    cap = max_thickness_ratio * W
    if cap <= W:
        return dict(ok=False, reason="thickness cap smaller than brick depth",
                    delta=0.0)
    if N < 2:
        return dict(ok=False, reason="need at least 2 courses", delta=0.0)
    hi = min(W * 0.999, (cap - W) / (N - 1))
    if max_overhang is not None:
        hi = min(hi, max_overhang)
    if hi <= 0.0:
        return dict(ok=False, reason="thickness cap too small for overhang",
                    delta=0.0)
    lo = 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if delta_feasible(mid, L, W, H, N, f_t, density, g, max_thickness_ratio,
                          max_overhang=max_overhang, bond_offset=bond_offset):
            lo = mid
        else:
            hi = mid
    return dict(ok=True, delta=lo)


def analyse(delta, L, W, H, N, f_t, density=BRICK_DENSITY, g=G,
            max_thickness_ratio=3.0, mortar_t=DEFAULT_MORTAR_T, t_m=None,
            bond_offset=0.0):
    """Full result set for a given uniform overhang."""
    t_m = mortar_t if t_m is None else t_m
    joints = joint_results(delta, L, W, H, N, density, g, bond_offset=bond_offset)
    max_stress = max(j["stress"] for j in joints)
    governing = max(joints, key=lambda j: j["stress"])
    toppling_margin = min(j["contact"] / 2.0 - j["ecc"] for j in joints)
    thickness = W + (N - 1) * delta
    centres = course_centres(delta, N)
    heights = [i * (H + t_m) + H / 2.0 for i in range(N)]
    angle = math.degrees(math.atan(delta / (H + t_m))) if (H + t_m) > 0 else 0.0
    return dict(
        delta=delta, joints=joints, max_stress=max_stress,
        governing_joint=governing["k"], stress_ratio=max_stress / f_t,
        toppling_margin=toppling_margin, thickness=thickness,
        thickness_cap=max_thickness_ratio * W,
        centres=centres, heights=heights, lean_angle=angle,
    )


def centres_from_deltas(deltas):
    """Overhang-direction centre (x) of each course from per-course overhangs.

    ``deltas[i]`` is the overhang of course i relative to course i-1
    (deltas[0] = 0).  Returns the list of course centres with course 0 at x=0.
    """
    x = [0.0]
    for d in deltas[1:]:
        x.append(x[-1] + d)
    return x


def profile_joint_results(xs, L, W, H, density=BRICK_DENSITY, g=G, bond_offset=0.0):
    """Evaluate every bed joint for an arbitrary (curved) course layout.

    ``xs`` are the course centres along the overhang direction.  Joint k carries
    the bricks at and above course k (k = 0 is the base joint under the first
    course, with full contact depth W).  Returns the same fields as
    ``joint_results`` (P, e, contact, moment, Z, net tensile stress).  ``bond_offset``
    reduces the effective joint width at every intermediate (running-bond) joint.
    """
    N = len(xs)
    Wb = brick_weight(L, W, H, density, g)
    w_int = max(L - bond_offset, 0.0)
    res = []
    contact = W
    ecc = sum(xs) / N
    P = N * Wb
    Z = L * contact * contact / 6.0
    res.append(
        dict(k=0, P=P, ecc=ecc, contact=contact,
             moment=P * ecc, Z=Z,
             stress=joint_stress(P, ecc, contact, L))
    )
    for k in range(1, N):
        n_above = N - k
        xbar = sum(xs[k:]) / n_above
        dk = xs[k] - xs[k - 1]
        contact = W - dk
        cc = (xs[k - 1] + xs[k]) / 2.0
        ecc = xbar - cc
        P = n_above * Wb
        Z = w_int * contact * contact / 6.0 if contact > 0 else 0.0
        res.append(
            dict(k=k, P=P, ecc=ecc, contact=contact,
                 moment=P * ecc, Z=Z,
                 stress=joint_stress(P, ecc, contact, w_int))
        )
    return res


def profile_feasible(deltas, L, W, H, f_t, density, g, cap, tol=1e-9,
                     max_overhang=None, bond_offset=0.0):
    """True if a curved corbel profile is stable (strength/toppling/cap)."""
    xs = centres_from_deltas(deltas)
    joints = profile_joint_results(xs, L, W, H, density, g, bond_offset=bond_offset)
    if max(j["stress"] for j in joints) > f_t:
        return False
    for j in joints:
        if j["ecc"] > j["contact"] / 2.0 + tol:
            return False
    if any(d >= W for d in deltas[1:]):
        return False
    if max_overhang is not None and any(d > max_overhang for d in deltas):
        return False
    if W + (xs[-1] - xs[0]) > cap + 1e-9:
        return False
    return True


def max_corbel_profile(L, W, H, N, f_t, density=BRICK_DENSITY, g=G,
                       max_thickness_ratio=3.0, max_overhang=None,
                       bond_offset=0.0, tol=1e-9):
    """Maximise the total lean with a curved (ogee) corbel profile.

    Starts from the block-stacking optimum (each course pushed so the centre of
    mass of everything above it sits exactly over the leading edge of the brick
    below), then scales it back only as far as the strength, per-brick overhang
    and thickness-cap limits require.  Returns a per-course overhang vector.
    """
    cap = max_thickness_ratio * W
    if cap <= W:
        return dict(ok=False, reason="thickness cap smaller than brick depth")
    if N < 2:
        return dict(ok=False, reason="need at least 2 courses")
    Wb = brick_weight(L, W, H, density, g)

    # Block-stacking optimal layout, built top-down.  The centre of mass of the
    # upper stack is placed at the leading edge of the course being added below.
    xs = [0.0] * N
    xs[N - 1] = 0.0
    weight = Wb
    cm = xs[N - 1]
    for k in range(N - 2, -1, -1):
        xs[k] = cm - W / 2.0
        w2 = weight + Wb
        cm = (weight * cm + Wb * xs[k]) / w2
        weight = w2
    x0 = xs[0]
    shape = [0.0] + [xs[i] - xs[i - 1] for i in range(1, N)]

    def scaled(s):
        return [0.0] + [s * d for d in shape[1:]]

    def feasible(s):
        return profile_feasible(scaled(s), L, W, H, f_t, density, g, cap, tol,
                                max_overhang=max_overhang, bond_offset=bond_offset)

    s = 1.0
    if not feasible(1.0):
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if feasible(mid):
                lo = mid
            else:
                hi = mid
        s = lo
    return dict(ok=True, deltas=scaled(s))


def analyse_profile(deltas, L, W, H, f_t, density=BRICK_DENSITY, g=G,
                    max_thickness_ratio=3.0, mortar_t=DEFAULT_MORTAR_T, t_m=None,
                    bond_offset=0.0):
    """Full result set for a curved corbel profile (per-course overhangs)."""
    t_m = mortar_t if t_m is None else t_m
    N = len(deltas)
    xs = centres_from_deltas(deltas)
    joints = profile_joint_results(xs, L, W, H, density, g, bond_offset=bond_offset)
    max_stress = max(j["stress"] for j in joints)
    governing = max(joints, key=lambda j: j["stress"])
    toppling_margin = min(j["contact"] / 2.0 - j["ecc"] for j in joints)
    thickness = W + (xs[-1] - xs[0])
    heights = [i * (H + t_m) + H / 2.0 for i in range(N)]
    dx = xs[-1] - xs[0]
    dz = heights[-1] - heights[0]
    angle = math.degrees(math.atan2(dx, dz)) if dz > 0 else 0.0
    return dict(
        deltas=deltas, centres=xs, joints=joints, max_stress=max_stress,
        governing_joint=governing["k"], stress_ratio=max_stress / f_t,
        toppling_margin=toppling_margin, thickness=thickness,
        thickness_cap=max_thickness_ratio * W,
        heights=heights, lean_angle=angle,
    )


if __name__ == "__main__":
    import Rhino.Geometry as rg

    MM = 0.001   # input length unit -> SI metres

    # -- GhPython inputs (lengths in mm) -----------------------------------
    brickL = globals().get("brickL", 190.0) * MM
    brickW = globals().get("brickW", 90.0) * MM
    brickH = globals().get("brickH", 57.0) * MM
    mortar_tensile = globals().get("mortar_tensile", 0.2e6)
    n_courses = globals().get("n_courses", 10)
    max_thickness_ratio = globals().get("max_thickness_ratio", 3.0)
    brick_density = globals().get("brick_density", BRICK_DENSITY)
    mortar_t = globals().get("mortar_t", 10.0) * MM

    n = int(n_courses)
    ratio = float(max_thickness_ratio)
    profile = globals().get("profile", 0)

    # Brick length (brickL) runs along the wall (the running/bond direction),
    # so the wall thickness / corbel overhang direction is the brick WIDTH.
    contact = brickW      # overhang (lean) depth  = brick width
    run_width = brickL    # wall-run (bond) width  = brick length

    # Default per-brick overhang cap = the overhang-direction brick width.
    mo = globals().get("max_overhang", None)
    max_overhang = (contact if mo is None else float(mo) * MM)

    # Running-bond offset of the odd courses = half the brick LENGTH, along the
    # wall run (Y).  The reduced bearing width is used in the stability check.
    bo = globals().get("bond_offset", None)
    bond_offset = (0.5 * run_width if bo is None else float(bo) * MM)

    if int(profile) == 1:
        res = max_corbel_profile(run_width, contact, brickH, n, mortar_tensile,
                                 density=brick_density,
                                 max_thickness_ratio=ratio,
                                 max_overhang=max_overhang,
                                 bond_offset=bond_offset)
    else:
        res = max_uniform_overhang(run_width, contact, brickH, n, mortar_tensile,
                                   density=brick_density,
                                   max_thickness_ratio=ratio,
                                   max_overhang=max_overhang,
                                   bond_offset=bond_offset)

    if not res["ok"]:
        Ok = False
        LeanAngle = 0.0
        Overhangs = [0.0] * n
        PerBrickX = [0.0] * n
        StressRatio = 1.0
        GoverningJoint = 0
        Frames = []
        _reason = res["reason"]
    else:
        if int(profile) == 1:
            analys = analyse_profile(res["deltas"], run_width, contact, brickH,
                                     f_t=mortar_tensile, density=brick_density,
                                     max_thickness_ratio=ratio, mortar_t=mortar_t,
                                     bond_offset=bond_offset)
            Overhangs = [d * M_TO_MM for d in analys["deltas"]]
        else:
            analys = analyse(res["delta"], run_width, contact, brickH, n,
                             mortar_tensile, density=brick_density,
                             max_thickness_ratio=ratio, mortar_t=mortar_t,
                             bond_offset=bond_offset)
            Overhangs = [0.0] + [analys["delta"] * M_TO_MM] * (n - 1)

        Ok = True
        LeanAngle = analys["lean_angle"]
        PerBrickX = [c * M_TO_MM for c in analys["centres"]]
        StressRatio = analys["stress_ratio"]
        GoverningJoint = analys["governing_joint"]

        Frames = []
        for i in range(n):
            y_off = bond_offset * M_TO_MM if i % 2 == 1 else 0.0
            org = rg.Point3d(analys["centres"][i] * M_TO_MM, y_off,
                             analys["heights"][i] * M_TO_MM)
            Frames.append(rg.Plane(org, rg.Vector3d(1.0, 0.0, 0.0),
                                   rg.Vector3d(0.0, 1.0, 0.0)))
