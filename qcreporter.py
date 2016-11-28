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

    def setStat(self, stat):
        self.stat = stat

    def setVersion(self, ver):
        self.version = ver

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
        self.outputSummary(io)
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
        io.write(".summary-table {padding:5px;border:1px solid #eeeeee;width:800px}"+"\n")
        io.write(".col1 {text-align:right;padding:5px;padding-right:20px;color:#666666;}"+"\n")
        io.write(".col2 {text-align:left;padding:5px;padding-left:20px;color:#332299;}"+"\n")
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
        idx = 1
        io.write("<li class='menu-item'><a href='#summary'>1, AfterQC summary</a> </li>\n")
        for figure in self.figures:
            idx += 1
            io.write("<li class='menu-item'><a href='#" + formatDivID(figure[0]) + "'>" + str(idx) + ", " + figure[0] + "</a> </li>\n")
        io.write("</ul></div>\n")

    def outputRow(self, io, k, v):
        io.write("<tr><td class='col1'>" + str(k) + "</td><td class='col2'>" + str(v) + "</td></tr>\n")

    def getSequencing(self):
        ret = str(self.stat["summary"]["readlen"])
        if self.stat["command"]["read2_file"] != None:
            ret = "2*" + ret + " pair end"
        else:
            ret += " single end"
        return ret

    def outputSummary(self, io):
        io.write("<div class='figure-div'>\n")
        io.write("<div class='figure-title'><a name='summary'>1, AfterQC summary</a></div>\n")
        io.write("<table class='summary-table'>\n")
        self.outputRow(io, "AfterQC Version:", self.version)
        self.outputRow(io, "sequencing:", self.getSequencing())
        if self.stat["command"]["read2_file"] != None:
            self.outputRow(io, "estimated seq error:", str(self.stat["overlap"]["error_rate"]*100) + "%")
        self.outputRow(io, "total reads:", self.stat["summary"]["total_reads"])
        self.outputRow(io, "filtered out reads:", str(self.stat["summary"]["bad_reads"]) + " <font color='#aaaaaa'>(" + str(100.0 * float(self.stat["summary"]["bad_reads"])/float(self.stat["summary"]["total_reads"])) + "%)</font>")
        self.outputRow(io, "total bases:", self.stat["summary"]["total_bases"])
        self.outputRow(io, "filtered out bases:", str(self.stat["summary"]["total_bases"] - self.stat["summary"]["good_bases"]) + " <font color='#aaaaaa'>(" + str(100.0 * float(self.stat["summary"]["total_bases"] - self.stat["summary"]["good_bases"])/float(self.stat["summary"]["total_bases"])) + "%)</font>")
        self.outputRow(io, "auto trimming", "front:" + str(self.stat["command"]["trim_front"]) + ", tail:" + str(self.stat["command"]["trim_tail"]) + " (use <font color='#aaaaaa'>-f0 -t0</font> to disable)")
        io.write("</table>\n")
        io.write("</div>\n")

    def outputFigures(self, io):
        io.write("<div id='figures'>")
        idx = 1
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
