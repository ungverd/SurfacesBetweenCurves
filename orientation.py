from typing import List, Literal
import bpy
import mathutils

def verify_selection(cur: bpy.types.Curve) -> Literal[False] | List[bpy.types.BezierSplinePoint]:
    points = [p for s in cur.splines for p in s.bezier_points if p.select_control_point]
    if len(points) == 2:
        if points[0].handle_left_type in ('ALIGNED', 'AUTO') and \
           points[0].handle_right_type in ('ALIGNED', 'AUTO') and \
           points[1].handle_left_type in ('ALIGNED', 'AUTO') and \
           points[1].handle_right_type in ('ALIGNED', 'AUTO'):
            return points
    elif len(points) == 1:
        if points[0].handle_left_type in ('VECTOR', 'FREE') or \
           points[0].handle_right_type in ('VECTOR', 'FREE'):
            return points
    raise ValueError(str(len(points)))
    return False

class CreateCoordsFromHandles(bpy.types.Operator):
    """My Creating Coords From Handles Script"""      # Use this as a tooltip for menu items and buttons.
    bl_idname = "curve.create_coords"        # Unique identifier for buttons and menu items to reference.
    bl_label = "Create coords"         # Display name in the interface.
    bl_options = {'REGISTER', 'UNDO'}  # Enable undo for the operator.

    def execute(self, context):        # execute() is called when running the operator.
        space = context.area.spaces[0]
        active = context.active_object
        cur = active.data
        points = verify_selection(cur)
        if not points:
            raise ValueError("error extractiong points")
        if len(points) == 2:
            vec1: mathutils.Vector = points[0].handle_left - points[0].co
            vec2: mathutils.Vector = points[1].handle_left - points[1].co
        else:
            vec1: mathutils.Vector = points[0].handle_left - points[0].co
            vec2: mathutils.Vector = points[0].handle_right - points[0].co
        vec3 = vec1.cross(vec2)
        vec3.normalize()
        rotation = mathutils.Vector((0,0,1)).rotation_difference(vec3)
        position = points[0].co + vec3
        view_matrix = mathutils.Matrix.LocRotScale(position, rotation, None)
        view_matrix.invert()
        space.region_3d.view_matrix = view_matrix
        bpy.ops.transform.create_orientation(name="handles", use=True)
        context.scene.transform_orientation_slots[0].custom_orientation.matrix = rotation.to_matrix()

        return {'FINISHED'}            # Lets Blender know the operator finished successfully.

'''class VIEW3D_MT_edit_curve_context_menu(bpy.types.Menu):
    bl_label = ''

    # Leave empty for compatibility.
    def draw(self, context):
        pass'''

'''def draw(self, context):
    layout = self.layout

    active = context.active_object
    cur = active.to_curve(context.evaluated_depsgraph_get())
    if verify_selection(cur):
        layout.separator()
        layout.operator(CreateCoordsFromHandles.bl_idname)

if __name__ == '__main__':
    # Register menu only if it doesn't already exist.
    rcmenu = getattr(bpy.types, "VIEW3D_MT_edit_curve_context_menu", None)
    if rcmenu is None:
        bpy.utils.register_class(VIEW3D_MT_edit_curve_context_menu)
        rcmenu = VIEW3D_MT_edit_curve_context_menu

    # Retrieve a python list for inserting draw functions.
    draw_funcs = rcmenu._dyn_ui_initialize()
    draw_funcs.append(draw)'''

def menu_func(self, context):
    self.layout.operator(CreateCoordsFromHandles.bl_idname)

def register():
    bpy.utils.register_class(CreateCoordsFromHandles)
    bpy.types.VIEW3D_MT_edit_curve_context_menu.append(menu_func)  # Adds the new operator to an existing menu.

def unregister():
    bpy.utils.unregister_class(CreateCoordsFromHandles)
if __name__ == '__main__':
    register()