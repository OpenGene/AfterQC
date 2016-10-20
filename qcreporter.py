import os,sys

def formatDivID(str):
    str = str.replace(" ", "-")
    str = str.replace(".", "-")
    str = str.replace("/", "-")
    return str

class QCReporter:

    def __init__(self):
        self.figures = []

    def addFigure(self, title, figure, div='', summary=''):
        self.figures.append((title, figure, div, summary))

    def output(self, filename):
        io = open(filename, "w")
        self.outputHeader(io)

    def outputHeader(self, io):
        io.write("<HTML>\n")
        io.write("<HEAD>\n")
        io.write('<script src="https://cdn.plot.ly/plotly-1.2.0.min.js"></script>\n')
        self.outputCSS(io)
        io.write("</HEAD>\n")
        io.write("<BODY>\n")
        io.write("<DIV id='container'>\n")
        self.outputMenu(io)
        self.outputFigures(io)
        io.write("</DIV>\n")
        self.outputPlotly(io)
        io.write("</BODY>\n")
        io.write("</HTML>")
        io.flush()

    def outputCSS(self, io):
        io.write("<style type=\"text/css\">"+"\n")
        io.write("#menu {text-align:left;}"+"\n")
        io.write(".menu-item{font-size:14px;padding:4px;}"+"\n")
        io.write("#container {text-align:center;padding-left:30px;}"+"\n")
        io.write(".figure-title {color:#bbbbbb;font-size:30px;padding:10px;text-align:left;}"+"\n")
        io.write(".figure-div {margin-top:40px;text-align:center}"+"\n")
        io.write(".plotly-div {width:800;height:600;text-align:center;}"+"\n")
        io.write("li {color:#666666;font-size:15px;border:0px;}"+"\n")
        io.write("</style>"+"\n")

    def outputPlotly(self, io):
        io.write("<script type=\"text/javascript\">"+"\n")
        for figure in self.figures:
            content = figure[1]
            div = figure[2]
            if div != '':
                io.write(content)
                io.write("\n")
        io.write("</script>"+"\n")

    def outputMenu(self, io):
        io.write("<div id='menu'><ul>\n")
        idx = 0
        for figure in self.figures:
            idx += 1
            io.write("<li class='menu-item'><a href='#" + formatDivID(figure[0]) + "'>" + str(idx) + ", " + figure[0] + "</a> </li>\n")
        io.write("</ul></div>\n")

    def outputFigures(self, io):
        io.write("<div id='figures'>")
        idx = 0
        for figure in self.figures:
            title = figure[0]
            content = figure[1]
            div = figure[2]
            summary = figure[3]
            idx += 1
            io.write("<div class='figure-div'>\n")
            io.write("<div class='figure-title'><a name='" + formatDivID(title) + "'>" + str(idx) + ", " + title + "</a></div>\n")
            if div == '':
                # an image figure
                io.write("<div class='figure'><img src='" + content + "'></div>\n")
            else:
                # a plotly figure
                io.write("<div id='" + div + "' class='plotly-div'></div>\n")
            io.write("<div class='figure-summary'>" + summary + "</div>\n")
            io.write("</div>\n")
        io.write("</div>\n")
