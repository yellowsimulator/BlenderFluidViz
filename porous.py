def main(context):
    """
    Adding a cube to the scene
    """
    bpy.ops.mesh.primitive_cube_add()
    
    

class PorousMedium(bpy.types.Operator):
    """
    """
    bl_idname = "myops.add_porous_medium"
    bl_label = "Add porous medium"


    def execute(self, context):
        main(context)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(PorousMedium)


def unregister():
    bpy.utils.unregister_class(PorousMedium)


if __name__=="__main__":
    register()
    bpy.ops.myops.add_porous_medium()