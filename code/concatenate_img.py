import cv2
import numpy as np
folder_path = "/home/zhangganlin/Desktop/Kombu/code/scenes/project/disney/metallic/"
img1 = cv2.imread(folder_path+"test_disney_metallic_0.png")
img2 = cv2.imread(folder_path+"test_disney_metallic_025.png")
img3 = cv2.imread(folder_path+"test_disney_metallic_050.png")
img4 = cv2.imread(folder_path+"test_disney_metallic_075.png")
img5 = cv2.imread(folder_path+"test_disney_metallic_100.png")


image = np.concatenate((img1,img2,img3,img4,img5),axis=1)

cv2.imwrite(folder_path+"concat_metallic.png",image)

# cv2.imshow('image',image)
# cv2.waitKey(0)