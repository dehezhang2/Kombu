import mitsuba as mi
mi.set_variant('scalar_rgb')
scene = mi.load_file("test_conductor_mitsuba.xml")
img = mi.render(scene)
mi.Bitmap(img).write('test_conductor_mitsuba.exr')

