#Video Vibrometry Matlab App

Based on [phase-based video vibrometry technique by MIT CSAIL](http://people.csail.mit.edu/nwadhwa/phase-video/).
This repo also acts as supplementary material for [my stereo vibrometry paper](https://www.mdpi.com/2075-4450/13/4/310/htm), where phase-based video vibrometry is used on multiple synchronized videos to extract 3D vibrational signals from black widow spiderwebs.

App that coordinates the use of phase-based optical flow to get pixel displacement estimates from video motion files.  

Instructions for use:
1) Open Matlab
2) `appdesigner` in Matlab command line
3) Open VideoVibrometry.mlapp
4) Run app
5) Load video in tab 1
6) Crop video to subsection to be analyzed
7) Choose sampling areas of video or choose to analyze entire video (or subregion if cropped in previous step) in tab 2
8) Process
9) Generate plots in tab 3
10) Save resulting displacement data
