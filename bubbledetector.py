 #!/usr/bin/env python
 
import os,sys
import fastq
import fileinput
import multiprocessing
from PIL import Image, ImageDraw
from circledetector import CircleDetector
    
########################### BubbleDetector
class BubbleDetector:
    
    def __init__(self, xmax, ymax, xmin, ymin, draw):
        #TODO: these can be set as parameters
        self.scope = 500.0
        self.cellSize = 300.0
        self.minPoly = 30
        self.minPointOfCluster = 150

        self.filename = ""
        self.rawRecords = []
        self.records = []
        self.meanNeighbour = 0.0
        self.gridX = 0
        self.gridY = 0
        self.grid = {}

        self.totalLabel = 0
        self.labels = []

        self.circles = []
        self.tile = 0
        self.lane = 0
        self.xmax = xmax
        self.ymax = ymax
        self.xmin = xmin
        self.ymin = ymin
        self.needDraw = draw
        self.initGrid()
    
    def detect(self):            
        #pass 1 filtering
        self.calcMeanCount()
        self.calcDensity()
        self.filterRecord(4)
        
        #pass2 filtering
        #self.calcMeanCount()
        #self.calcDensity()
        #self.filterRecord(100)

        #merge all records for clustering
        self.mergeRecord()
        
        #clustering by region grow
        #give every record a label
        self.clustering()
        self.filterCluster()
        
        self.detectCircles()
        
        if len(self.circles)>0:
            print("raw records: " + str(len(self.rawRecords)))
            print("filtered records: " + str(len(self.records)))
            print("raw labels: " + str(self.totalLabel))
            print("filtered labels:" + str(len(self.labels)))
            print("bubbles:" + str(len(self.circles)))
            print(self.circles)
            if self.needDraw:
                print("drawing image...")
                self.draw(self.records)
                
        return self.circles
            
    def initGrid(self):
        self.gridX = int(self.xmax/self.cellSize)+1
        self.gridY = int(self.ymax/self.cellSize)+1
        
        for gx in xrange(self.gridX):
            self.grid[gx]={}
            for gy in xrange(self.gridY):
                self.grid[gx][gy] = []
                
    def calcMeanCount(self):
        totalCount = 0
        for gx in xrange(self.gridX):
            for gy in xrange(self.gridY):
                for rec in self.grid[gx][gy]:
                    totalCount += rec[3]
        percent = (self.scope*self.scope*3.1415926)/float(self.xmax*self.ymax)
        self.meanNeighbour = percent * totalCount
                
    def countNeighbour(self, rec):
        
        x = rec[0]
        y = rec[1]
        
        gxmin = int((x-self.scope)/self.cellSize)
        gxmin = max(0,gxmin)
        gxmax = int((x+self.scope)/self.cellSize) + 1
        gxmax = min(self.gridX - 1, gxmax)
        
        gymin = int((y-self.scope)/self.cellSize)
        gymin = max(0,gymin)
        gymax = int((y+self.scope)/self.cellSize) + 1
        gymax = min(self.gridY - 1, gymax)
                
        for ngx in xrange(gxmin, gxmax+1):
            for ngy in xrange(gymin, gymax+1):
                for nrec in self.grid[ngx][ngy]:
                    nx = nrec[0]
                    ny = nrec[1]
                    if (nx - x)*(nx - x) + (ny - y)*(ny - y) < self.scope * self.scope:
                        count = nrec[3]
                        rec[4] += count
                            
    def calcDensity(self):
        for gx in xrange(self.gridX):
            for gy in xrange(self.gridY):
                for rec in self.grid[gx][gy]:
                    rec[4] = 0
                    self.countNeighbour(rec)
            
    def mergeRecord(self):
        self.records = []
        for gx in xrange(self.gridX):
            for gy in xrange(self.gridY):
                self.records = self.records + self.grid[gx][gy]
                
    def isNeighbour(self, rec1, rec2):
        x1 = rec1[0]
        y1 = rec1[1]
        x2 = rec2[0]
        y2 = rec2[1]
        if (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) < self.scope * self.scope:
            return True
        else:
            return False
            
    def filterCluster(self):
        labels = xrange(1, self.totalLabel+1)
        numbers = {}
        for label in labels:numbers[label] = 0
        
        for rec in self.records:
            label = rec[6]
            numbers[label] += 1
        
        self.labels = []
        for label in labels:
            if numbers[label] >= self.minPointOfCluster and label not in self.labels:
                self.labels.append(label)
         
        for i in range(len(self.records))[::-1]:
            rec = self.records[i]
            if rec[6] not in self.labels:
                self.records.pop(i)
        
                
    def clustering(self):
        queue = []
    
        while True:
            curX = 0
            curY = 0
            
            #find the first record with no label
            findNoLabel = False
            for rec in self.records:
                label = rec[6]
                if label == 0:
                    findNoLabel = True
                    self.totalLabel += 1
                    rec[6] = self.totalLabel
                    queue.append(rec)
                    break
             
            #done       
            if findNoLabel == False:
                break
                
            while len(queue) != 0:
                #label all neighbours, and put them to queue
                seed = queue.pop()
                curX = seed[0]
                curY = seed[1]
                curLabel = seed[6]
                for rec in self.records:
                    label = rec[6]
                    if label == 0:
                        if self.isNeighbour(rec, seed):
                            rec[6] = curLabel
                            queue.append(rec)
                            
        #sort records by label
        self.records.sort(key=lambda x: (x[6]))
                    
    
    def filterRecord(self, fold):
        minNeighbour = self.meanNeighbour * fold
        
        for gx in xrange(self.gridX):
            for gy in xrange(self.gridY):
                for i in range(len(self.grid[gx][gy]))[::-1]:
                    rec = self.grid[gx][gy][i]
                    neighbour = rec[4]
                    if neighbour < minNeighbour:
                        self.grid[gx][gy].pop(i)
    
    def setFilename(self, filename):
        self.filename  = filename
    
    def setTile(self, tile):
        self.tile  = tile 
        
    def setLane(self, lane):
        self.lane  = lane 
    
    def loadRecords(self, records):
        for r in records:
            surface = r[1]
            x = r[5]
            y = r[6]
            base = r[7]
            count = r[8]
               
            #place holders
            density = 0.0
            label = 0
            neighbour = 0
                
            #we only care polyT or polyG here
            if base == ord('G') and count >= self.minPoly:
                rec = [x, y, surface, count, neighbour, density, label]
                gx = int(x/self.cellSize)
                gy = int(y/self.cellSize)
                self.rawRecords.append(rec)
                self.grid[gx][gy].append(rec)
                    
    def loadRecordsFromFile(self, filename):
        self.filename = filename
        with open(filename) as f:
            rows = f.readlines()
            for row in rows[1:]:
                r = row.split(",")
                surface = int(r[1])
                x = int(r[5])
                y = int(r[6])
                base = int(r[7])
                count = int(r[8])
                
                #place holders
                density = 0.0
                label = 0
                neighbour = 0
                
                #we only care polyT or polyG here
                if base == ord('G') and count >= self.minPoly:
                    rec = [x, y, surface, count, neighbour, density, label]
                    gx = int(x/self.cellSize)
                    gy = int(y/self.cellSize)
                    self.rawRecords.append(rec)
                    self.grid[gx][gy].append(rec)
    
    def detectCircles(self):
        lastLabel = -1
        labelRecords = []
        
        for rec in self.records:
            label = rec[6]
            
            if lastLabel != label and lastLabel != -1:
                cd = CircleDetector(labelRecords, self.xmax, self.ymax, self.xmin, self.ymin)
                labelCircles = cd.detect()
                self.circles = self.circles + labelCircles
                #clear it
                labelRecords = []
                
            labelRecords.append(rec)
            lastLabel = label
            
        if len(labelRecords) != 0:
            cd = CircleDetector(labelRecords, self.xmax, self.ymax, self.xmin, self.ymin)
            labelCircles = cd.detect()
            self.circles = self.circles + labelCircles
        
        #write lane and tile    
        for c in self.circles:
            c[4] = self.lane
            c[5] = self.tile
        
    def draw(self, polyRecords):
            
        #define the colors of A T C G
        colors = {}
        for c in xrange(100):
            red = (c * 29791)%128 + 127
            green = (c * 67571)%128 + 127
            blue = (c * 87571)%128 + 127
            colors[c] = [red, green, blue]

        xMax = self.xmax
        yMax = self.ymax
        countMax = 50
        
        #TODO: make this a parameter
        tileImageScale = 0.1
        tileImageWidth = int(tileImageScale*xMax) + 4
        tileImageHeight = int(tileImageScale*yMax) + 4
        
        #draw pixels        
        tileImageData = [[0,0,0] for x in xrange(tileImageWidth * tileImageHeight)]
        for r in polyRecords:
            x = r[0]
            y = r[1]
            surface = r[2]
            count = r[3]
            label = r[6]
            
            inCircle = False
            for circle in self.circles:
                if circle[3] == label:
                    inCircle = True
                    break
                    
            if inCircle == False:
                continue
                            
            #calc the alpha
            alpha = 1.0
            blendColor = colors[label]
            
            #########################################
            #calc the tile image pixel pos
            tileImagePixelX = int(x * tileImageScale)
            tileImagePixelY = int(y * tileImageScale)
            
            #get original tile image pixel data
            pixel = tileImageData[tileImageWidth * tileImagePixelY + tileImagePixelX]
            
            #blend
            for c in xrange(3):
                pixel[c] = blendColor[c]
                pixel[c] = min(255, pixel[c])
                
            #write back the tile image  pixel
            tileImageData[tileImageWidth * tileImagePixelY + tileImagePixelX] = pixel
        
        
        img =  Image.new("RGB", (tileImageWidth, tileImageHeight), "black")
        img.putdata([tuple(x) for x in tileImageData])
        
        draw = ImageDraw.Draw(img)
         #draw circles
        for circle in self.circles:
            centerX = circle[0] * tileImageScale
            centerY = circle[1] * tileImageScale
            radius = circle[2] * tileImageScale
            draw.ellipse((centerX-radius, centerY-radius, centerX+radius, centerY+radius))
            
        img.save(self.filename + "." + str(len(self.circles)) + ".png")

#test
def main():        
    dashplace = 0
    for argv in sys.argv[1:]:
        print(argv)
        xmax = 0
        ymax = 0
        xmin = 999999999
        ymin = 999999999
        with open(argv) as f:
            rows = f.readlines()
            for row in rows[1:]:
                r = row.split(",")
                x = int(r[5])
                y = int(r[6])
                xmin = min(xmin,x)
                ymin = min(ymin,y)
                xmax = max(xmax, x)
                ymax = max(ymax, y)
        bd = BubbleDetector(xmax, ymax, xmin,ymin, True)
        bd.loadRecordsFromFile(argv)
        bd.detect()
