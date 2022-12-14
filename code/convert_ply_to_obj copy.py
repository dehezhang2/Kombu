import os
in_folder = "/home/zhangganlin/Desktop/Kombu/main_scene-20221124T111915Z-001/main_scene/meshes"
out_folder = "/home/zhangganlin/Desktop/Kombu/main_scene-20221124T111915Z-001/main_scene/obj_meshes"
names_ply = os.listdir(in_folder)
names = [ply[:-4]for ply in names_ply]
names_obj = [out_folder+"/"+name+".obj" for name in names]
names_ply = [in_folder+"/"+name+".ply"for  name in names]

for i in range(len(names)):
    os.system("ctmconv \""+names_ply[i]+"\" \""+names_obj[i]+"\"")