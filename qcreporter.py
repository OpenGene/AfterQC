import os,sys

def formatDivID(str):
    str = str.replace(" ", "-")
    str = str.replace(".", "-")
    str = str.replace("/", "-")
    return str

class QCReporter:

    def __init__(self):
        self.figures = []

    def addFigure(self, title, figure, summary=''):
        self.figures.append((title, figure, summary))

    def output(self, filename):
        io = open(filename, "w")
        self.outputHeader(io)

    def outputHeader(self, io):
        io.write("<HTML>")
        io.write("<HEAD>")
        self.outputCSS(io)
        io.write("</HEAD>")
        io.write("<BODY>")
        io.write("<DIV id='container'>")
        self.outputMenu(io)
        self.outputFigures(io)
        io.write("</DIV>")
        io.write("</BODY>")
        io.write("</HTML>")

    def outputCSS(self, io):
        io.write("<style type=\"text/css\">"+"\n")
        io.write("#menu {text-align:left;}"+"\n")
        io.write(".menu-item{font-size:20px;padding:5px;}"+"\n")
        io.write("#container {text-align:center;}"+"\n")
        io.write(".figure-title {color:#ff6600;font-weight:bold;font-size:20px;padding:10px;}"+"\n")
        io.write(".figure-div {margin-top:40px;}"+"\n")
        io.write("li {color:#666666;font-size:15px;border:0px;}"+"\n")
        io.write("</style>"+"\n")

    def outputMenu(self, io):
        io.write("<div id='menu'><ul>")
        for figure in self.figures:
            io.write("<li class='menu-item'><a href='#" + formatDivID(figure[0]) + "'>" + figure[0] + "</a> </li>")
        io.write("</ul></div>")

    def outputFigures(self, io):
        io.write("<div id='figures'>")
        for figure in self.figures:
            io.write("<div class='figure-div'>")
            io.write("<div class='figure-title'><a name='" + formatDivID(figure[0]) + "'>" + figure[0] + "</a></div>")
            io.write("<div class='figure'><img src='" + figure[1] + "''></div>")
            io.write("<div class='figure-summary'>" + figure[2] + "</div>")
            io.write("</div>")
        io.write("</div>")
