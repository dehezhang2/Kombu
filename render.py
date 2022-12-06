import mitsuba as mi
mi.set_variant('scalar_rgb')
scene = mi.load_file("earth/earth/earth.xml")
img = mi.render(scene)
mi.Bitmap(img).write('earth.exr')

