#This project contains the following 3 .cpp files

#1.minray.cpp - The original code which produces the image
#2.minraycopy.cpp - For changing input size
#3.minraycopy2.cpp - For changing number of threads

#Commands to compile and execute on Windows:-

g++ minray.cpp -o minray -fopenmp -std=gnu++11
minray.exe > minray.ppm