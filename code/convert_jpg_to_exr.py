import os
in_folder = "/home/zhangganlin/Desktop/Kombu/main_scene-20221124T111915Z-001/main_scene/textures"
out_folder = "/home/zhangganlin/Desktop/Kombu/main_scene-20221124T111915Z-001/nori/texture"
names_ply = os.listdir(in_folder)
names = [ply[:-4]for ply in names_ply]
names_obj = [out_folder+"/"+name+".exr" for name in names]
names_ply = [in_folder+"/"+name+".jpg"for  name in names]

for i in range(len(names)):
    os.system("convert \""+names_ply[i]+"\" \""+names_obj[i]+"\"")