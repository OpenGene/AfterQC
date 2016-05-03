 #!/usr/bin/env python
 
import os,sys
import time
import fastq
import fileinput
import multiprocessing
import math
    
########################### circleDetector
class CircleDetector:
    
    def __init__(self, allRecords, xmax, ymax, xmin, ymin):
        self.records = allRecords
        self.xmax = xmax
        self.ymax = ymax
        self.xmin = xmin
        self.ymin = ymin
        self.records = []
        self.radiusThreshold = 800
        self.areaRatioThreshold = 0.85
    
    def detect(self):    
        circles = []
        if self.isInCorner():
            return circles
    
        c = self.detectDirectly()
        
        boundArea = self.boundRectArea()
        circleArea = self.circleAreaInTile(c, self.xmax, self.ymax)
        
        areaRatio = float(boundArea) / circleArea
        if areaRatio < self.areaRatioThreshold:
            return circles
        
        if c[2]>self.radiusThreshold:
            circles.append(c)
        return circles
        
    def circleAreaInTile(self, c, width, height):
        #we need to deal with the cases that circle in on edge(s) of tile
        x = c[0]
        y = c[1]
        r = c[2]
        circleArea = math.pi * r * r
        
        xout = 0
        yout = 0
        if x-r-self.xmin < 0:xout = -(x-r-self.xmin)
        if y-r-self.ymin < 0:yout = -(y-r-self.ymin)
        if x+r-width > xout :xout = x+r-width
        if y+r-height > yout :yout = y+r-height

        #considering a circle with r=1.0
        if xout>0 or yout>0:
            out = xout + yout
            outProp = out/r
            #this means center is outside the tile, in this case we turn to calc the inner area
            if out>r :
                outProp = 2.0 - outProp
            #calc the angle oriented from the center
            angle = math.acos(1.0 - outProp)
            fanArea = angle/2.0
            triangleArea = math.sin(angle) * (1.0 - outProp) * 0.5
            #up and down, so double it
            outareaProp = (fanArea - triangleArea) * 2.0
            if out>r :
                outareaProp = math.pi - outareaProp
            return circleArea * (1.0 - outareaProp/math.pi)

        return circleArea


    def boundRectArea(self):
        left = self.xmax
        top = self.ymax
        right = 0
        bottom = 0
        
        for rec in self.records:
            x = rec[0]
            y = rec[1]
            left = min(left, x)
            top = min(top, y)
            right = max(right, x)
            bottom = max(bottom, y)
            
        return (right - left) * (bottom - top)
        
    def isInCorner(self):
        
        cornerCount = 0
        thresholdX = 0.1 * self.xmax
        thresholdY = 0.1 * self.ymax
        for rec in self.records:
            x = rec[0]
            y = rec[1]
            if (x < thresholdX or x > self.xmax - thresholdX) and (y < thresholdY or y > self.ymax - thresholdY):
                cornerCount += 1
        
        if float(cornerCount) / len(self.records) > 0.1:
            return True
        else:
            return False
        
    def detectOnEdge(self):
        
        topnumber = 500
        top = [[0,0,0] for x in xrange(topnumber)]
        centerX = 0
        centerY = 0
        
        #find the pair of two points with max distance
        #set its center as circle center
        for i in xrange(len(self.records)):
            rec1 = self.records[i]
            x1 = rec1[0]
            y1 = rec1[1]
            for j in xrange(i):
                rec2 = self.records[j]
                x2 = rec2[0]
                y2 = rec2[1]
                dist = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
                pos = topnumber
                for t in xrange(topnumber):
                    if dist < top[t][2]:
                        pos = t+1
                        break
                if pos > 0:
                    for m in xrange(1, pos):
                        top[m-1] = top[m]
                    cx = (x1+x2)/2.0
                    cy = (y1+y2)/2.0
                    top[pos-1] = [cx, cy, dist]
                    
        totalX = 0
        totalY = 0
        totalNum = 0
        for point in top:
            if top[2]>0:
                totalX = totalX + point[0]
                totalY = totalY + point[1]
                totalNum += 1
                    
        centerX  = totalX / float(totalNum)
        centerY  = totalY / float(totalNum)
        
        #calc radius
        maxDistance = 0.0
        for rec in self.records:
            x = rec[0]
            y = rec[1]
            distance = (centerX - x) * (centerX - x) + (centerY - y) * (centerY - y)
            if distance > maxDistance:
                maxDistance = distance
                
        radius = math.sqrt(float(maxDistance))
        label = self.records[0][6]
        
        #tile and lane is a place holder
        lane = 0
        tile = 0
        circle = [centerX, centerY, radius, label, lane, tile]
        return circle
    
    #old method, deprecated
    def detectDirectly(self):
        totalX = 0
        totalY = 0
        
        #calc center
        for rec in self.records:
            x = rec[0]
            y = rec[1]
            totalX = totalX + x
            totalY = totalY + y
        
        recNum = len(self.records)  
        centerX = totalX/float(recNum)
        centerY = totalY/float(recNum)
        
        #calc radius
        maxDistance = 0
        for rec in self.records:
            x = rec[0]
            y = rec[1]
            distance = (centerX - x) * (centerX - x) + (centerY - y) * (centerY - y)
            if distance > maxDistance:
                maxDistance = distance
                
        r = math.sqrt(maxDistance)
        x = centerX
        y = centerY
        
        if x-r-self.xmin < 0 or y-r-self.ymin < 0 or x+r-self.xmax > 0 or y+r-self.ymax > 0:
            if len(self.records)<5000:
                return self.detectOnEdge()
        
        label = self.records[0][6]
        
        #tile and lane is a place holder
        lane = 0
        tile = 0
        circle = [centerX, centerY, r, label, lane, tile]
        return circle

def main():
    cd = CircleDetector([], 10.0, 10.0)
    c = [5, -0.5, 1.0]
    print(cd.circleAreaInTile(c, 10.0, 10.0))
