"""Repeats an optimised corbelled brick column along the wall length.

Takes the column frames produced by ``brick_corbel_column`` and tiles a copy of
the whole column for each column position along the wall length.  The column is
one brick long (the brick length runs along the wall, the frame Y axis), so
columns are spaced ``brickL + joint_width`` apart along the wall run.

Inputs (all lengths in mm, matching the model unit; geometry is output in mm):
    frames:        Column centre frames (planes) from brick_corbel_column
    brickL:        Brick length - runs along the wall (frame Y direction) [mm]
    brickW:        Brick width - the wall thickness / overhang direction [mm]
    brickH:        Brick height [mm]
    wall_length:   Straight-wall run length [mm]; NOT used for a curved wall.
                   Only applies when no base_curve is supplied.  For a curve the
                   wall length is the measured ARC length of the base_curve
    joint_width:   Mortar joint gap [mm] (used only when bricks_per_course = 0)
    bricks_per_course: Exact number of bricks in every course (int, 0 = auto).
                   When > 0 the bricks are spread evenly along the curve so the
                   spacing AUTO-varies to fill it; then wall_length/joint_width
                   are not used for placement
    brick_gap:     Changeable horizontal gap between bricks [mm] (default =
                   joint_width); used only when bricks_per_course = 0
    min_gap:       Minimum gap between neighbouring bricks [mm] (default 3)
    max_lean:      Cap on the total overhang of the top brick from its column
                   base [mm] (default = half the brick length).  Columns whose
                   top rows would overhang more are flattened to this limit
    bond_offset:   Running-bond offset of the odd courses along the wall run
                   [mm]; overrides the column's value.  Lower it (e.g. 45mm) so
                   the combined lean + bond overhang stays under half a brick
                   length (default = the column's bond offset)
    base_curve:    Rhino curve from curve_generate.py (a sine wave) that defines
                   course 0: one column is placed for each sample along the curve,
                   and each column is rotated so its run direction follows the
                   curve tangent so the bricks stay in contact.  The overhang
                   points toward the sine centre and tapers to 0 at the
                   zero-crossings so the wall stays connected through direction
                   changes.  If bricks would collide, the gap is widened
                   automatically until there is no overlap.  If unset, a straight
                   wall is tiled from wall_length instead
Outputs:
    WallFrames:    Centre frame of every brick in the wall
    WallBricks:    Box (Rhino.Geometry.Box) for every brick in the wall
    NColumns:      Number of columns placed along the wall length
    WallLength:    Actual wall length [mm] - the base_curve arc length for a
                   curved wall, or the straight-wall run when no curve is used
"""

__author__ = "aa9663"
__version__ = "2023.09.11"

import math
import Rhino.Geometry as rg

X = rg.Vector3d(1.0, 0.0, 0.0)
Y = rg.Vector3d(0.0, 1.0, 0.0)


def flatten(data):
    """Flatten a DataTree, a [DataTree], or nested lists into one flat list."""
    if hasattr(data, "Branches"):
        items = []
        for b in data.Branches:
            items.extend(list(b))
        return items
    if hasattr(data, "AllData"):
        return list(data.AllData())
    if isinstance(data, (list, tuple)):
        out = []
        for it in data:
            out.extend(flatten(it))
        return out
    return [data]


def _origin(frame):
    """Origin of a frame, or the point itself if a Point3d is passed."""
    return frame.Origin if hasattr(frame, "Origin") else frame


def _axes(frame):
    """Return (xaxis, yaxis); default to world axes for a plain point."""
    if hasattr(frame, "XAxis"):
        return frame.XAxis, frame.YAxis
    return X, Y


def _rect_poly(cx, cy, ux, uy, hl, hw):
    """Convex 2D polygon (plan footprint) of a brick centred at (cx, cy).

    ``u = (ux, uy)`` is the unit run direction (half-length hl), the normal is
    taken perpendicular (half-width hw).
    """
    vx, vy = -uy, ux
    return [
        (cx + ux * hl + vx * hw, cy + uy * hl + vy * hw),
        (cx + ux * hl - vx * hw, cy + uy * hl - vy * hw),
        (cx - ux * hl - vx * hw, cy - uy * hl - vy * hw),
        (cx - ux * hl + vx * hw, cy - uy * hl + vy * hw),
    ]


def _poly_area(poly):
    area = 0.0
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return abs(area) * 0.5


def _clip_poly(subject, clip):
    """Sutherland-Hodgman clip of a convex subject by a convex clip polygon."""
    # Normalise the "clip" polygon to counter-clockwise so the inside test works.
    sa = 0.0
    for i in range(len(clip)):
        x1, y1 = clip[i]
        x2, y2 = clip[(i + 1) % len(clip)]
        sa += x1 * y2 - x2 * y1
    if sa < 0:
        clip = clip[::-1]

    def _inside(p, a, b):
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) >= -1e-12

    def _inter(p1, p2, a, b):
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = a
        x4, y4 = b
        den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if abs(den) < 1e-12:
            return p1
        tt = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
        return (x1 + tt * (x2 - x1), y1 + tt * (y2 - y1))

    output = subject
    n = len(clip)
    for i in range(n):
        a = clip[i]
        b = clip[(i + 1) % n]
        inp = output
        output = []
        if not inp:
            break
        s = inp[-1]
        for e in inp:
            if _inside(e, a, b):
                if not _inside(s, a, b):
                    output.append(_inter(s, e, a, b))
                output.append(e)
            elif _inside(s, a, b):
                output.append(_inter(s, e, a, b))
            s = e
    return output


def _overlap_area(poly_a, poly_b):
    """Plan overlapping area (mm^2) of two convex brick footprints."""
    clip = _clip_poly(poly_a, poly_b)
    return _poly_area(clip) if len(clip) >= 3 else 0.0


def _chop_curve(curve, n=256):
    """Sample a curve into (params, cumulative_arc_length, total_length)."""
    dom = curve.Domain
    params = [dom.T0 + (dom.T1 - dom.T0) * i / n for i in range(n + 1)]
    cum = [0.0]
    for i in range(1, n + 1):
        p1 = curve.PointAt(params[i - 1])
        p2 = curve.PointAt(params[i])
        dd = math.sqrt((p2.X - p1.X) ** 2 + (p2.Y - p1.Y) ** 2
                       + (p2.Z - p1.Z) ** 2)
        cum.append(cum[-1] + dd)
    return params, cum, cum[-1]


def _param_at_arc(params, cum, s):
    """Clamp s to [0, total] and return the curve parameter at arc length s."""
    total = cum[-1]
    if total <= 0:
        return params[0]
    if s <= 0:
        return params[0]
    if s >= total:
        return params[-1]
    lo, hi = 0, len(cum) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if cum[mid] < s:
            lo = mid + 1
        else:
            hi = mid
    if lo == 0:
        return params[0]
    seg = (s - cum[lo - 1]) / max(cum[lo] - cum[lo - 1], 1e-12)
    return params[lo - 1] + seg * (params[lo] - params[lo - 1])


if __name__ == "__main__":
    frames = flatten(globals().get("frames", []))
    brickL = globals().get("brickL", 190.0)
    brickW = globals().get("brickW", 90.0)
    brickH = globals().get("brickH", 57.0)
    wall_length = globals().get("wall_length", 2000.0)
    joint_width = globals().get("joint_width", 10.0)
    bricks_per_course = globals().get("bricks_per_course", 0)
    brick_gap = globals().get("brick_gap", None)
    min_gap = globals().get("min_gap", 3.0)
    max_lean = globals().get("max_lean", None)
    bond_offset = globals().get("bond_offset", None)     # running-bond offset [mm]
    base_curve = globals().get("base_curve", None)

    # Cap on the total overhang (top brick's offset from its column base) so the
    # top rows don't stick out more than this (default = half the brick length).
    if max_lean is None:
        max_lean = 0.5 * brickL

    step = brickL + joint_width

    # Coerce base_curve to a single Rhino curve if one is supplied.
    curve = None
    if base_curve is not None:
        if isinstance(base_curve, (list, tuple)) and base_curve:
            base_curve = base_curve[0]
        try:
            import rhinoscriptsyntax as rs
            curve = rs.coercecurve(base_curve)
        except Exception:
            curve = None
        if curve is None:
            curve = base_curve

    if curve is not None and step > 0:
        # Course 0 follows the input curve: lay a column for each sample point,
        # rotated so the run direction (brick length) follows the curve tangent
        # and the overhang/lean follows the curve normal (outward).
        length = float(curve.GetLength())
        dom = curve.Domain

        # Centre of the sine (the line it oscillates around) = mid-point of the
        # X range, sampled along the curve.
        xs = []
        for i in range(65):
            p = curve.PointAt(dom.T0 + (dom.T1 - dom.T0) * i / 64.0)
            xs.append(p.X)
        center_x = 0.5 * (min(xs) + max(xs))
        amp = 0.5 * (max(xs) - min(xs))

        xr = rg.Interval(-brickW / 2.0, brickW / 2.0)
        yr = rg.Interval(-brickL / 2.0, brickL / 2.0)
        zr = rg.Interval(-brickH / 2.0, brickH / 2.0)

        # Arc-length table for placing columns at a given distance along the curve.
        params_arc, cum_arc, total_arc = _chop_curve(curve)
        n_courses = len(frames)

        # Gap / count control.  If bricks_per_course > 0 it is used exactly;
        # otherwise the count is derived from the curve length and the gap.
        gap = (brick_gap if brick_gap is not None else joint_width)
        ncols_req = int(bricks_per_course) if bricks_per_course else 0

        # Maximum total overhang of the column profile (top-course offset from its
        # base), in mm - used to flatten the top-heavy rows.
        max_oX = 0.0
        for fr in frames:
            max_oX = max(max_oX, _origin(fr).X)

        def col_data(s):
            """Build one column at arc length s into (pt, tan, normal, lean_f,
            origins, course_footprint_polys)."""
            t = _param_at_arc(params_arc, cum_arc, s)
            pt = curve.PointAt(t)
            tan = rg.Vector3d(curve.TangentAt(t))
            tan.Z = 0.0
            if tan.Length < 1e-9:
                tan = rg.Vector3d(0.0, 1.0, 0.0)
            tan.Unitize()
            normal = rg.Vector3d(-tan.Y, tan.X, 0.0)
            if (pt.X - center_x) * normal.X > 0:      # overhang toward the centre
                normal = rg.Vector3d(-normal.X, -normal.Y, 0.0)
            if amp < 1e-9:
                lean_f = 1.0
            else:
                lean_f = min(max(abs(pt.X - center_x) / amp, 0.0), 1.0)
            # Flatten the corbel: cap the total overhang so the top rows do not
            # stick out more than `max_lean`.
            if max_oX > 1e-9:
                lean_f = min(lean_f, max_lean / max_oX)
            origins = []
            polys = []
            for fr in frames:
                o = _origin(fr)
                lean_off = o.X * lean_f
                # Running-bond offset: use the column's value unless overridden.
                oY = o.Y
                if bond_offset is not None and oY != 0.0:
                    oY = float(bond_offset)
                org = rg.Point3d(
                    pt.X + normal.X * lean_off + tan.X * oY,
                    pt.Y + normal.Y * lean_off + tan.Y * oY,
                    o.Z,
                )
                origins.append(org)
                polys.append(_rect_poly(org.X, org.Y, tan.X, tan.Y,
                                        brickL / 2.0, brickW / 2.0))
            return pt, tan, normal, lean_f, origins, polys

        def clears(s_a, s_b, tol=1e-6):
            """True if the columns at s_a and s_b keep at least ``min_gap``
            between every brick (all courses) - i.e. the bricks never overlap
            even when grown by min_gap/2 each side."""
            da = col_data(s_a)
            db = col_data(s_b)
            hl = brickL / 2.0 + min_gap / 2.0
            hw = brickW / 2.0 + min_gap / 2.0
            for k in range(n_courses):
                ea = _rect_poly(da[4][k].X, da[4][k].Y, da[1].X, da[1].Y, hl, hw)
                eb = _rect_poly(db[4][k].X, db[4][k].Y, db[1].X, db[1].Y, hl, hw)
                if _overlap_area(ea, eb) > tol:
                    return False
            return True

        def any_collision(s_list, tol=1e-6):
            """True if ANY two columns (adjacent or not) overlap in plan on ANY
            course - the final brute-force check for truly collision-free output.
            """
            datas = [col_data(si) for si in s_list]
            for i in range(len(datas)):
                da = datas[i]
                for j in range(i + 1, len(datas)):
                    db = datas[j]
                    for k in range(n_courses):
                        ea = _rect_poly(da[4][k].X, da[4][k].Y, da[1].X, da[1].Y,
                                        brickL / 2.0, brickW / 2.0)
                        eb = _rect_poly(db[4][k].X, db[4][k].Y, db[1].X, db[1].Y,
                                        brickL / 2.0, brickW / 2.0)
                        if _overlap_area(ea, eb) > tol:
                            return True
            return False

        def placement_clear(s_list, tol=1e-6):
            """True if the placement keeps min_gap between every adjacent pair
            AND has no overlap between ANY two columns on ANY course."""
            for i in range(1, len(s_list)):
                if not clears(s_list[i - 1], s_list[i], tol):
                    return False
            return not any_collision(s_list, tol)

        WallFrames = []
        WallBricks = []
        NColumns = 0
        WallLength = length
        if frames:
            if ncols_req > 0:
                # Fixed brick count: each column is placed at the smallest arc
                # length that keeps at least `min_gap` from the previous one, so
                # the only gaps that grow are between colliding bricks.  If a
                # brick genuinely cannot fit (curve too short for the count + gap)
                # the count is reduced so the wall never collides.
                for trial in range(ncols_req):
                    ncols = ncols_req - trial
                    if ncols < 2:
                        NColumns = 0
                        break
                    s = [0.0] * ncols
                    s[0] = 0.0
                    fits = True
                    for i in range(1, ncols):
                        lo = s[i - 1]
                        hi = total_arc
                        if not clears(s[i - 1], hi):
                            fits = False
                            break                    # reduce the count & retry
                        for _ in range(45):
                            mid = 0.5 * (lo + hi)
                            if clears(s[i - 1], mid):
                                hi = mid
                            else:
                                lo = mid
                        s[i] = hi
                    if not fits:
                        continue
                    # Reject if any adjacent pair is under min_gap or any two
                    # columns overlap on any course.
                    if not placement_clear(s):
                        continue
                    # Centre the wall only if the shifted placement is still clear
                    # (shifting along the curve changes the lean/normal geometry).
                    span = s[-1] - s[0]
                    if span <= total_arc:
                        shift = 0.5 * (total_arc - span) - s[0]
                        cand = [x + shift for x in s]
                        if placement_clear(cand):
                            s = cand
                    for si in s:
                        pt, tan, normal, lean_f1, origins, polys = col_data(si)
                        for org in origins:
                            plane = rg.Plane(org, normal, tan)
                            WallFrames.append(plane)
                            WallBricks.append(rg.Box(plane, xr, yr, zr))
                    NColumns = ncols
                    break
            else:
                # Auto mode: derive the count from the curve length and gap, place
                # the columns at a constant spacing, widening the gap on collision.
                for _attempt in range(25):
                    spacing = brickL + gap
                    ncols = int(length // spacing) if spacing > 0 else 0
                    if ncols < 2:
                        ncols = 0
                        break
                    start_arc = 0.5 * (total_arc - (ncols - 1) * spacing)
                    if start_arc < 0:
                        start_arc = 0.0

                    frames_out = []
                    bricks_out = []
                    origin_xs = []
                    tangents = []
                    for c in range(ncols):
                        s = start_arc + c * spacing
                        pt, tan, normal, lean_f1, origins, polys = col_data(s)
                        origin_xs.append(origins)
                        tangents.append((tan.X, tan.Y))
                        for org in origins:
                            plane = rg.Plane(org, normal, tan)
                            frames_out.append(plane)
                            bricks_out.append(rg.Box(plane, xr, yr, zr))

                    max_overlap = 0.0
                    hl = brickL / 2.0 + min_gap / 2.0
                    hw = brickW / 2.0 + min_gap / 2.0
                    for k in range(n_courses):
                        for c in range(ncols - 1):
                            ea = _rect_poly(origin_xs[c][k].X, origin_xs[c][k].Y,
                                            tangents[c][0], tangents[c][1], hl, hw)
                            eb = _rect_poly(origin_xs[c + 1][k].X, origin_xs[c + 1][k].Y,
                                            tangents[c + 1][0], tangents[c + 1][1], hl, hw)
                            ov = _overlap_area(ea, eb)
                            if ov > max_overlap:
                                max_overlap = ov

                    WallFrames = frames_out
                    WallBricks = bricks_out
                    NColumns = ncols
                    if max_overlap <= 1e-6:
                        break
                    gap += max(1.0, max_overlap / max(brickW, 1e-9))
                    if gap > total_arc:
                        break
    else:
        ncols = int(wall_length // step) if step > 0 else 0
        WallFrames = []
        WallBricks = []
        WallLength = ncols * step if ncols > 0 else 0.0
        if frames and ncols > 0:
            # Brick length runs along the wall run (frame Y); the brick width is
            # the wall thickness / overhang direction (frame X).
            xr = rg.Interval(-brickW / 2.0, brickW / 2.0)
            yr = rg.Interval(-brickL / 2.0, brickL / 2.0)
            zr = rg.Interval(-brickH / 2.0, brickH / 2.0)

            # Centre the wall around x=0 (across the required length)
            total = ncols * step
            y0 = -total / 2.0

            for c in range(ncols):
                yoff = y0 + c * step + step / 2.0
                for fr in frames:
                    org = _origin(fr)
                    xaxis, yaxis = _axes(fr)
                    org = rg.Point3d(org) + rg.Vector3d(0.0, yoff, 0.0)
                    plane = rg.Plane(org, xaxis, yaxis)
                    WallFrames.append(plane)
                    WallBricks.append(rg.Box(plane, xr, yr, zr))
        NColumns = ncols
