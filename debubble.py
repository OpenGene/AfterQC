 #!/usr/bin/env python
 
import os,sys
from optparse import OptionParser
import time
import fastq
from bubbleprocesser import BubbleProcesser

def parseCommand():
    usage = "usage: %prog <input_files> [options]"
    version = "%prog 1.1"  
    parser = OptionParser(usage = usage, version = version) 
    parser.add_option("-p", "--poly_max", dest = "poly_max", default = 20, type = "int",
        help = "if exists one polyX(polyG means GGGGGGGGG...), and its length is >= poly_max, then this read is considered probabily in a bubble")
    parser.add_option("-o", "--output", dest = "output", default = "bubble_detect",
        help = "folder to store the csv and image files")
    parser.add_option("-i", "--input", dest = "input", default = ".",
        help = "folder storing input fastq files")
    parser.add_option("-d", "--draw", dest = "draw", default = False,
        help = "specify whether draw the pictures or not")
    return parser.parse_args()
    
def writeCircles(circles, outdir):
    outfile = open(os.path.join(outdir, "circles.csv"), "w")
    outfile.write("x,y,radius,lane,tile\n")
    for c in circles:
        outfile.write(str(c[0])+",")
        outfile.write(str(c[1])+",")
        outfile.write(str(c[2])+",")
        outfile.write(str(c[4])+",")
        outfile.write(str(c[5])+"\n")
    
def oldmain():
    time1 = time.time()
    (options, args) = parseCommand()
        
    dashplace = 0
    for argv in sys.argv[1:]:
        dashplace += 1
        if argv.startswith("-"):
            break
    if dashplace < 1:
        print("usage: %prog <input_files> [options]")
        print("use -h option to see help")
        sys.exit(1)
        
    bp = BubbleProcesser(sys.argv[1:dashplace], options.poly_max, options.output, options.draw)
    bp.run()
    
    time2 = time.time()
    print('Time used: ' + str(time2-time1))
    
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
    print("detected bubbles:")
    print(circles)
    writeCircles(circles, output)
    
def main():
    time1 = time.time()
    (options, args) = parseCommand()
        
    debubbleDir(options.input, options.poly_max, options.output, options.draw)
    
    time2 = time.time()
    print('Time used: ' + str(time2-time1))

if __name__  == "__main__":
    main()