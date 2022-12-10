import mitsuba as mi
mi.set_variant('scalar_rgb')
scene = mi.load_file("test_bumpmap_mitsuba.xml")
img = mi.render(scene)
mi.Bitmap(img).write('test_bumpmap_mitsuba.exr')

