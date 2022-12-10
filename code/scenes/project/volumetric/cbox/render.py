import mitsuba as mi
mi.set_variant('scalar_rgb')
scene = mi.load_file("cbox_vol_path_mats_misuba.xml")
img = mi.render(scene)
mi.Bitmap(img).write('cbox_vol_path_mats_misuba.exr')

