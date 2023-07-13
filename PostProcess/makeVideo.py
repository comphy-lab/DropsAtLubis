# Python script intended to convert a set of .png images to a video

import cv2
import os
# Variable components (Change imageFolder and videoName for each case)
# Case numenclature: 'videoCase#ChangedParameterValue'
# Example: 'videoCase1000Ohw1e-3.mp4' 
# See ListOfSimulations for more infromation in parant folder. 
imageFolder = '/media/vatsal/Spreading/TestPrecursorFilm/3cases/1001/Video'
videoName = 'videoCase1000Ohw1e-3.mp4'

# Sort images in the folder by .png (video will be in same folder)
images = sorted([img for img in os.listdir(imageFolder) if img.endswith('.png')])
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 

# Create path for video 
videoPath = os.path.join(imageFolder, videoName)
videoWriter = cv2.VideoWriter(videoPath, fourcc, 30.0, (1920, 1080))

# for loop which iterates each image and adds it to the video
for imgName in images:
    imgPath = os.path.join(imageFolder, imgName)
    img = cv2.imread(imgPath)
    # if condition to reshape the image to fit video format
    if img is None:
        print(f'Error loading image: {imgPath}')
        continue
    if img.shape != (1920, 1080):
        img = cv2.resize(img, (1920, 1080))
    for i in range(1): 
        videoWriter.write(img)
videoWriter.release()

# Checks if the video path exist otherwise outputs confirmation
if not os.path.exists(videoPath):
    print(f'Error creating video file: {videoPath}')
else:
    print(f'Video file created: {videoPath}')
