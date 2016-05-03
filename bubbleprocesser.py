 #!/usr/bin/env python
 
import os,sys
import re
import time
import fastq
import fileinput
import multiprocessing
from PIL import Image, ImageDraw
from multiprocessing import Process, Queue
from bubbledetector import BubbleDetector
    
########################### BubbleProcesser
class BubbleProcesser:

    __input = None
    __output = None
    __poly_x_num = None
    __drawImage = False
    xmax = 0
    ymax = 0
    xmin = 0
    ymin = 0
    
    def __init__(self, input, poly_x_num, output, drawImage):
        self.__input = input
        self.__output = output
        self.__poly_x_num = poly_x_num
        self.__drawImage = drawImage
        self.__polyRecords = []
        self.circles = []
        
        
    def run(self):
        
        if len(self.__input) == 0:
            print("No fastq files")
            return self.circles
        
        #create queues to store the result
        queues = {}
        for fname in self.__input:
            queues[fname]  = Queue()
        
        #create a job for each file
        jobs = [Process(target = self.statFile, args = (fname,queues[fname])) for fname in self.__input]
        for job in jobs: job.start()
        
        #wait for the result
        results = [queues[fname].get() for fname in self.__input]
        for job in jobs: job.join()
        
        #merge the result
        for r in results:
            self.__polyRecords = self.__polyRecords + r
            
        print("finished polyX stat for all files")
        
        self.calcMaxMin()
        
        #sort by surface
        self.__polyRecords.sort(key=lambda x: x[1])
        
        #sort by tileID
        self.__polyRecords.sort(key=lambda x: (x[9]%10000))
        
        #sort by lane
        self.__polyRecords.sort(key=lambda x: x[0])
        
        print("write records to poly_X.csv")
        self.writeToFile()
        
        print("process records by tile")
        self.processByTile()
        
        if(self.__drawImage):
            print("draw images")
            self.drawImages()
            
        return self.circles
        
    def calcMaxMin(self):
        maxValue = [0 for x in xrange(10)]
        for r in self.__polyRecords:
            for i in xrange(10): 
                maxValue[i] = max(maxValue[i], r[i])
        
        minValue = [maxValue[x] for x in xrange(10)]
        for r in self.__polyRecords:
            for i in xrange(10): 
                minValue[i] = min(minValue[i], r[i])        
        
        self.xmax = maxValue[5]
        self.ymax = maxValue[6]
        self.xmin = minValue[5]
        self.ymin = minValue[6]
        
    def drawImages(self):
        
        if not os.path.exists(self.__output):
            os.makedirs(options.output)
            
        if not os.path.exists(os.path.join(self.__output, "image_by_camera")):
            os.makedirs(os.path.join(self.__output, "image_by_camera"))
        
        if len(self.__polyRecords) == 0:
            return
            
        maxValue = [0 for x in xrange(10)]
        
        for r in self.__polyRecords:
            for i in xrange(10): 
                maxValue[i] = max(maxValue[i], r[i])
                    
        laneMin = self.__polyRecords[0][0]
        laneMax = self.__polyRecords[len(self.__polyRecords)-1][0]
        
        #separate records by lane
        laneRecords = {}
        for lane in xrange(laneMin, laneMax+1):
            laneRecords[lane]=[]
            
        for r in self.__polyRecords:
            lane = r[0]
            laneRecords[lane].append(r)
            
        #start jobs by lane
        jobs = [Process(target = self.draw, args = (laneRecords[laneRec], maxValue)) for laneRec in laneRecords]
        for job in jobs: job.start()
        
        #wait for the result
        for job in jobs: job.join()
            
    def detectBubbleForTile(self, tileRecords, tile_no, lane):
        if not os.path.exists(os.path.join(self.__output, "image_by_tile")):
            os.makedirs(os.path.join(self.__output, "image_by_tile"))
        bd = BubbleDetector(self.xmax, self.ymax, self.xmin, self.ymin, self.__drawImage)
        bd.loadRecords(tileRecords)
        bd.setLane(lane)
        bd.setTile(tile_no)
        bd.setFilename(os.path.join(self.__output, "image_by_tile", str(lane) + "_" + str(tile_no)))
        c = bd.detect()
        self.circles = self.circles + c
        
    def processByTile(self):
        if not os.path.exists(os.path.join(self.__output, "record_by_tile")):
            os.makedirs(os.path.join(self.__output, "record_by_tile"))
        
        tileRecords = []
        lastTileWithBothSurface = -1
        laneOfLastTile = -1
        
        for r in self.__polyRecords:
            lane = r[0]
            tile_no = r[9] % 10000
            
            #save this tile
            if tile_no != lastTileWithBothSurface:
                if lastTileWithBothSurface != -1:
                    self.writeRecords(os.path.join(self.__output, "record_by_tile", str(laneOfLastTile) + "_x" + str(lastTileWithBothSurface)+".csv"), tileRecords)
                    self.detectBubbleForTile(tileRecords, lastTileWithBothSurface, laneOfLastTile)
                    tileRecords = []
                lastTileWithBothSurface = tile_no
            
            if lane != laneOfLastTile:
                laneOfLastTile = lane    
            
            tileRecords.append(r)
        
        #the last tile is not written, so write it here    
        if lastTileWithBothSurface != -1:
            #save last tile
            self.writeRecords(os.path.join(self.__output, "record_by_tile", str(laneOfLastTile) + "_x" + str(lastTileWithBothSurface)+".csv"), tileRecords)
            self.detectBubbleForTile(tileRecords, lastTileWithBothSurface, laneOfLastTile)
        
    def draw(self, polyRecords, maxValue):
            
        #define the colors of A T C G
        colors = {}
        colors[ord('A')] = [255, 0, 0]
        colors[ord('T')] = [0, 255, 0]
        colors[ord('C')] = [0, 0, 255]
        colors[ord('G')] = [150,60,240]
        
        lane = polyRecords[0][0]
                
        laneMax = maxValue[0]
        surfaceMax = maxValue[1]
        swathMax = maxValue[2]
        cameraMax = maxValue[3]
        tileMax = maxValue[4]
        xMax = max(25920, maxValue[5])
        yMax = max(19440, maxValue[6])
        countMax = min(100, maxValue[8])
        
        #TODO: make these parameters
        cameraImageScale = 0.01
        tileGap = 0.1
        margin = 50
        
        cameraImageWidth = (swathMax * xMax * (1.0 + tileGap)) * cameraImageScale + 2*margin
        cameraImageHeight = (tileMax * yMax * (1.0 + tileGap)) * cameraImageScale + 2*margin
        cameraImageWidth = int(cameraImageWidth)
        cameraImageHeight = int(cameraImageHeight)
        
        #createimage data buffers
        cameraImageData = {}
        cameraImageDataCount = {}
        for camera in xrange(0, cameraMax+1):
            cameraImageData[camera] = [[0,0,0] for x in xrange(cameraImageWidth * cameraImageHeight)]
            cameraImageDataCount[camera] = 0
        
        #draw pixels        
        
        #remember that the records are sorted by tile, so if we meet a new tile, that means last tile is finished
        lastTile = -1
        for r in polyRecords:
            lane = r[0]
            surface = r[1]
            swath = r[2]
            camera = r[3]
            tile = r[4]
            x = r[5]
            y = r[6]
            base = r[7]
            count = r[8]
            tile_no = r[9] % 10000
            
            cameraImageDataCount[camera] += 1
                
            #calc the alpha
            alpha = float(count)/float(countMax)
            blendColor = colors[base]
            
            #########################################
            #calc the camera image pixel pos
            cameraImagePixelX = ((swath - 1) * (1.0 + tileGap) * xMax + x) * cameraImageScale + margin
            cameraImagePixelY = ((tile - 1) * (1.0 + tileGap) * yMax + y) * cameraImageScale + margin
            cameraImagePixelX = int(cameraImagePixelX)
            cameraImagePixelY = int(cameraImagePixelY)
            
            #calc the camera image pixel offset
            cameraImageOffset = cameraImagePixelY * cameraImageWidth + cameraImagePixelX
                
            #get original camera image pixel data
            pixel = cameraImageData[camera][cameraImageOffset]
            
            #blend
            for c in xrange(3):
                pixel[c] = int(alpha * blendColor[c] + (1.0 - alpha) * pixel[c])
                pixel[c] = min(255, pixel[c])
                
            #write back the camera image  pixel
            cameraImageData[camera][cameraImageOffset] = pixel
        
        #########################################
            
        #write camera images
        print("write images by each camera")
        for camera in xrange(0, cameraMax+1):
            #skip the tiles that has no data, which actually don't exist
            if cameraImageDataCount[camera] == 0:
                continue
            img =  Image.new("RGB", (cameraImageWidth, cameraImageHeight), "black")
            img.putdata([tuple(x) for x in cameraImageData[camera]])
            #draw circles
            for circle in self.circles:
                circleLane = circle[4]
                circleTile = str(circle[5])
                swath = int(circleTile[0])
                circleCamera = int(circleTile[1])
                tile = int(circleTile[2:])
                if circleLane != lane or camera != circleCamera:
                    continue
                x = circle[0]
                y = circle[1]
                lineWidth = 7
                radius = circle[2] * cameraImageScale - lineWidth/2
                centerX =  ((swath - 1) * (1.0 + tileGap) * xMax + x) * cameraImageScale + margin
                centerY = ((tile - 1) * (1.0 + tileGap) * yMax + y) * cameraImageScale + margin
                draw = ImageDraw.Draw(img)
                #draw concentrical circles to make the circle thicker
                for round in xrange(lineWidth):
                    draw.ellipse((centerX-radius, centerY-radius, centerX+radius, centerY+radius))
                    radius += 1.0
            img.save(os.path.join(self.__output, "image_by_camera", str(lane)+"_"+str(camera)+".png"))
            print("finished drawing lane " + str(lane) + " camera: " + str(camera))
            
    def writeToFile(self):
        if not os.path.exists(self.__output):
            os.makedirs(self.__output)
        self.writeRecords(os.path.join(self.__output, "poly_X.csv"), self.__polyRecords)
        
    def writeRecords(self, filename, records):
        outfile = open(filename, "w")
        outfile.write("lane,surface,swath,camera,tile,xpos,ypos,base,count, tile_no\n")
        if len(records) == 0:
            outfile.close()
            return
        for record in records:
            outfile.write(",".join([str(x) for x in record]) + "\n") 
        outfile.close()
    
    def statFile(self, filename, queue):
        if filename.endswith(".bam") or filename.endswith(".cram"):
            return self.statFileBam(filename, queue)
        else:
            return self.statFileFastq(filename, queue)
    
    def statFileBam(self, filename, queue):
        print("start: " + filename + "\n")
        #we apply lazy import here for pysam, since this module need to be installed manually and it is not used for most case
        import pysam
        reader = pysam.AlignmentFile(filename)
        records = []
        while True:
            try:
                read = reader.next()
            except StopIteration:
                break
            if read == None:
                break
            poly, count = self.countPoly(read.seq)
            
            if poly != None:
                
                #illumina sequence name line format
                #@<instrument>:<run number>:<flowcell ID>:<lane>:<tile_no>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
                
                items = read.qname.split(":")
                if len(items) < 7:
                    continue
                lane = items[3]
                tile_no = items[4]
                surface = tile_no[0]
                swath = tile_no[1]
                camera = tile_no[2]
                tile = tile_no[3:]
                xpos = items[5]
                ypos = items[6]
                record = [int(lane), int(surface), int(swath), int(camera), int(tile), int(xpos), int(ypos), ord(poly), count, int(tile_no)]
                records.append(record)

        queue.put(records)
        reader.close()
        print("finished " +filename + " with " + str(len(records)) + " polyX records")
                   
    def statFileFastq(self, filename, queue):
        print("start: " + filename + "\n")
        reader = fastq.Reader(filename)
        records = []
        pattern = re.compile(r'\S+\:\d+\:\S+\:\d+\:\d+\:\d+\:\d+')
        while True:
            read = reader.nextRead()
            if read == None:
                break
            poly, count = self.countPoly(read[1])
            
            if poly != None:
                
                #illumina sequence name line format
                #@<instrument>:<run number>:<flowcell ID>:<lane>:<tile_no>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
                
                match = pattern.search(read[0]);
                if not match:
                    continue

                items = match.group().split(":")
                if len(items) < 7:
                    continue
                lane = items[3]
                tile_no = items[4]
                surface = tile_no[0]
                swath = tile_no[1]
                camera = tile_no[2]
                tile = tile_no[3:]
                xpos = items[5]
                ypos = items[6]
                record = [int(lane), int(surface), int(swath), int(camera), int(tile), int(xpos), int(ypos), ord(poly), count, int(tile_no)]
                records.append(record)

        queue.put(records)
        print("finished " +filename + " with " + str(len(records)) + " polyX records")
        
    def countPoly(self, read):
        polyArray = ["A", "T", "C", "G"]
        for p in polyArray:
            if p*self.__poly_x_num in read:
                count = self.__poly_x_num
                pos = read.find(p*self.__poly_x_num)
                for c in read[pos + self.__poly_x_num:]:
                    if c == p:
                        count += 1
                    else:
                        break
                return p, count
        return None, 0