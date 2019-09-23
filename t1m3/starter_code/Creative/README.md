Submission
==========
A video recording of your creative scene should be placed in this directory
with one of these extensions:
 - .avi
 - .mov
 - .mpg
 - .mp4
 - .mkv

it should be named like this: `uni1234_t1m3.avi`.

For example, I'd submit as `jda2158_t1m3.avi`.

You should also include the XML scene file that generated your scene, which
should be named `uni1234_t1m3.xml`.

Be sure to use several of the new features we added this week! Here are a few
of them, to remind you:
 - Linearized implicit Euler
 - Vortex Forces

Creating Your Video
===================
After you've run cmake at least once, if you go into your build/ directory and
run ccmake .., you can hit the arrow keys to move to the option `USE_PNG`, and
then hit enter to toggle it from OFF to ON. Then, you can hit c to configure
your build to use this option. After that, recompiling with make should set up
your FOSSSim to create .png files when it runs.

It will name these frames, for exmaple ./pngs/frame00001.png. You therefore need
to be sure that a pngs/ folder exists in whatever directory you're in at the
time you run FOSSSim.

Every time you advance the simulation a frame, it will output a new frame file.

Once you're satisfied you've created all the frames you need, you can use ffmpeg
to generate a .avi file by following this guide:

	https://trac.ffmpeg.org/wiki/Create%20a%20video%20slideshow%20from%20images

We look forward to seeing what you'll make!
