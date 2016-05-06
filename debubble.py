 #!/usr/bin/env python
 
import os,sys
from optparse import OptionParser
import time
import fastq
from bubbleprocesser import BubbleProcesser
from util import *

def parseCommand():
    usage = "usage: %prog <input_files> [options]"
    version = "%prog 1.1"  
    parser = OptionParser(usage = usage, version = version) 
    parser.add_option("-p", "--poly_max", dest = "poly_max", default = 20, type = "int",
        help = "Min polyX to draw on tile images. Default is 20.")
    parser.add_option("-o", "--output", dest = "output", default = "bubble",
        help = "folder to store the csv and image files. Default is bubble.")
    parser.add_option("-i", "--input", dest = "input", default = ".",
        help = "folder storing input fastq files. Default is current dir.")
    parser.add_option("-d", "--draw", dest = "draw", default = "on",
        help = "specify whether draw the pictures or not. Default is on.")
    return parser.parse_args()
    
def writeCircles(circles, outdir):
    if not os.path.exists(outdir):
        return
    outfile = open(os.path.join(outdir, "circles.csv"), "w")
    outfile.write("x,y,radius,lane,tile\n")
    for c in circles:
        outfile.write(str(c[0])+",")
        outfile.write(str(c[1])+",")
        outfile.write(str(c[2])+",")
        outfile.write(str(c[4])+",")
        outfile.write(str(c[5])+"\n")
    
def debubbleDir(folder, poly_max, output, drawImage):
    files = os.listdir(folder)
    fastqs = []
    for f in files:
        path = os.path.join(folder, f)
        if os.path.isdir(path):
            continue
        if fastq.isFastq(path):
            fastqs.append(path)
            
    bp = BubbleProcesser(fastqs, poly_max, output, drawImage)
    circles = bp.run()
    if len(circles) >0 :
        print("detected bubbles:")
        print(circles)
    else:
        print("no bubble detected")
    writeCircles(circles, output)
    
def runInFolder():
    time1 = time.time()
    (options, args) = parseCommand()
    options.draw = parseBool(options.draw)
        
    debubbleDir(options.input, options.poly_max, options.output, options.draw)
    
    time2 = time.time()
    print('Time used in folder: ' + str(time2-time1))

def main():
    time1 = time.time()
    (options, args) = parseCommand()
        
    dashplace = 1
    for argv in sys.argv[1:]:
        if argv.startswith("-"):
            break
        dashplace += 1
    if dashplace <= 1:
        runInFolder()
        return
        
    bp = BubbleProcesser(sys.argv[1:dashplace], options.poly_max, options.output, options.draw)
    circles = bp.run()
    if len(circles) >0 :
        print("detected bubbles:")
        print(circles)
    else:
        print("no bubble detected")
    writeCircles(circles, options.output)
    
    time2 = time.time()
    print('Time used: ' + str(time2-time1))
