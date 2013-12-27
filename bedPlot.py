#!/usr/bin/env python
from Tkinter import *
import sys

# ---- classes -----
class CoordMapper:
   def __init__(self,zoom): self.zoom=zoom

# ---- callbacks ----

# ---- Main class ---
class SequenceFigure:
    def __init__(self,root, seqs, fontSize):
        self.root = root
        self.seqs = seqs
        self.fontSize = fontSize
        self.boxWidth=fontSize
        self.minx = 10
        self.miny = 10
        self.tickFont = ("Helvetica", str(fontSize), "")
        self.seqFont =  ("Courier", str(fontSize), "")
        sizeX, sizeY = 3000, 800
        self.setupCanvas(root, sizeX, sizeY)

    def setupCanvas(self, root, sizeX, sizeY):
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)
        canvas = Canvas(root, bg="yellow", width=800, height=400, scrollregion=(0,0,sizeX,sizeY))
        scrolly = Scrollbar(root, orient=VERTICAL, command=canvas.yview)
        scrolly.grid(row=0, column=1, sticky=N+S)
        scrollx = Scrollbar(root, orient=HORIZONTAL, command=canvas.xview)
        scrollx.grid(row=1, column=0, sticky=E+W)
        canvas.grid(row=0, column=0,sticky=N+S+E+W)
        canvas["xscrollcommand"] = scrollx.set
        canvas["yscrollcommand"] = scrolly.set

        self.canvas = canvas
        self.scrollx = scrollx
        self.scrolly = scrolly

    def zoomInc(self):
        self.canvas.delete(self.canvas.find_all())

    def zoomDec(self):
        self.canvas.delete(items=0)

    def drawRuler(self, miny, len, tickDist, labelDist):
        """ tickDist = distance between ticks (boxes) ,
            len = length of rulers in boxWidths, 
            labelDist = how many ticks between all labels
        """
        tickHeight = 5
        tickToLabelMargin = 5
        self.canvas.create_line(self.minx,miny,self.minx+len*self.boxWidth,miny)
        i = 0
        # draw ticks and labels
        #for x in range(self.minx,self.minx+len*self.boxWidth,tickDist*self.boxWidth):
        for c in range(0,len,tickDist):
            x = self.minx+c*self.boxWidth
            self.canvas.create_line(x,miny,x,miny+tickHeight)
            if labelDist != 0 and i % labelDist == 0:
                self.canvas.create_text(x,miny+tickHeight+tickToLabelMargin, text=str(i*tickDist), font=self.tickFont)
            i+=1

        return miny+tickHeight+self.fontSize

    def drawSeqs(self, lasty, seqs):
        for seqNo in range(0, len(seqs)):
            seq = seqs[seqNo]
            i = 0
            for x in range(self.minx, self.boxWidth*len(seq), self.boxWidth):
                self.canvas.create_text(x, lasty+seqNo, text=seq[i], font=self.seqFont)
                i+=1
            lasty+=self.fontSize
        lasty = lasty + len(seqs) * self.fontSize + self.fontSize
        return lasty
        
    def draw(self):
        lens = [len(s) for s in self.seqs]
        maxSeqLen = max(lens)
        trackMargin = 5

        lasty = self.drawRuler(self.miny, maxSeqLen, 10, 1)

        lasty += trackMargin
        lasty = self.drawSeqs (lasty, self.seqs)

class Actions:
    def __init__(self,seqFig,scrollx, scrolly, canvas):
        self.scrollx = scrollx
        self.scrolly = scrolly
        self.seqFig=seqFig
        self.canvas = canvas

    def zoomInc(self, event):
        self.seqFig.zoomInc()

    def zoomDec(self, event):
        self.seqFig.zoomDec()

    def quit(self, event):
        sys.exit(0)

    def scrollLeft(self, event):
        self.canvas.xview("scroll","-1","units")
    def scrollRight(self, event):
        self.canvas.xview("scroll","1","units")
    def jumpLeft(self, event):
        self.canvas.xview("scroll","-1","page")
    def jumpRight(self, event):
        self.canvas.xview("scroll","1","page")

# --- setup frame ----
def setupFrame(root, seqs):
    seqfig = SequenceFigure(root, seqs, 15)
    actions = Actions(seqfig, seqfig.scrollx, seqfig.scrolly, seqfig.canvas)

    # wire Keyboard shotcuts to actions
    root.bind("q", actions.quit)
    root.bind("<Key-Left>", actions.scrollLeft)
    root.bind("<Key-Right>", actions.scrollRight)
    root.bind("<Shift-Left>", actions.jumpLeft)
    root.bind("<Shift-Right>", actions.jumpRight)
    root.bind("+", actions.zoomInc)
    root.bind("-", actions.zoomDec)

    seqfig.draw()
    #canvas.pack(expand=1, side=LEFT,anchor=NW)
    #drawFigure(canvas, "", seqs, 15) 

# --- MAIN -----
def main():
    root = Tk()
    seqs = [10 * "aaacgtactgcatc", 10 * "aaacgtactgcatc", 10 * "aaccaccttgg"]
    setupFrame(root, seqs)
    root.mainloop()
    #drawing_area.bind("<Motion>", motion)
    #drawing_area.bind("<ButtonPress-1>", b1down)
    #drawing_area.bind("<ButtonRelease-1>", b1up)

def b1down(event):
    global b1
    b1 = "down"           # you only want to draw when the button is down
                          # because "Motion" events happen -all the time-

def b1up(event):
    global b1, xold, yold
    b1 = "up"
    xold = None           # reset the line when you let go of the button
    yold = None

def motion(event):
    if b1 == "down":
        global xold, yold
        if xold != None and yold != None:
            event.widget.create_line(xold,yold,event.x,event.y,smooth=TRUE)
                          # here's where you draw it. smooth. neat.
        xold = event.x
        yold = event.y

if __name__ == "__main__":
    main()
