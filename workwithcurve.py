import bpy
import mathutils
cur = bpy.data.curves['BezierCurve.001']
splines = cur.splines
vec = mathutils.Vector((0.01, 0.0, 0.0))
for s in splines:
    for p in s.bezier_points:
        if p.select_control_point:
            p.handle_left = p.co + vec
            p.handle_right = p.co - vec