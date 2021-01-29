class SVariant:
    def __init__(self, tool, line=None, chrom=None, pos=None, id=None, ref=None, end=None, gt=None, svlen=None, svtype=None, cipos1=None, cipos2=None, ciend1=None, ciend2=None, algorithms=None):
        self.tool = tool
        if(line is not None):
            self.parse_line(line)
        else:
            self.chrom = chrom
            self.pos = pos
            self.end = end
            self.gt = gt
            self.id = id
            self.ref = ref
            self.svlen = svlen
            self.svtype = svtype
            self.cipos1 = cipos1
            self.cipos2 = cipos2
            self.ciend1 = ciend1
            self.ciend2 = ciend2
            self.used = False
            self.algorithms = algorithms
    def printVcfLine(self):
        sv_len = self.svlen
        end = self.end
        if(self.svtype == "INS" or self.svtype == "DUP"):
            sv_len = abs(sv_len)
            if(self.svtype == "INS"):
                end = self.pos
        return '\t'.join((self.chrom, str(self.pos), self.id, self.ref, "<"+self.svtype+">",
                 ".", "PASS", "END="+str(end)+";SVLEN="+str(sv_len)+";SVTYPE="+self.svtype+";ALGORITHM="+self.algorithms+";CIPOS="+str(self.cipos1)+","+str(self.cipos2)+";CIEND="+str(self.ciend1)+","+str(self.ciend2), "GT", self.gt, "\n"))
    def parse_line(self, line):
        values = line.split("\t")
        self.chrom = values[0]
        self.pos = int(values[1])
        info = values[7].split(";")
        self.gt = values[9]
        self.ref = values[3]
        
        self.end = int(info[0].split("=")[1])
        self.svlen = info[1].split("=")[1]

        if(self.svlen == "."):
            self.svlen = self.end-self.pos
        self.svlen = int(self.svlen)

        self.svtype = self.parse_type(info[2].split("=")[1])
        if(self.svtype == "INS"):
            self.end = self.pos+abs(self.svlen)
        cipos = info[3].split("=")[1]
        ciend = info[4].split("=")[1]

        if(cipos == "."):
            cipos = "-10,10" # maybe other values? 0s?
        if(ciend == "."):
            ciend = "-10,10" # maybe other values? 0s?

        cipos = cipos.split(",")
        ciend = ciend.split(",")

        self.cipos1 = int(cipos[0])
        self.cipos2 = int(cipos[1])

        self.ciend1 = int(ciend[0])
        self.ciend2 = int(ciend[1])

        self.used = False
    def parse_type(self, svtype):
        if "del" in svtype.casefold():
            return "DEL"
        if "inv" in svtype.casefold():
            return "INV"
        if "ins" in svtype.casefold():
            return "INS"
        if "dup" in svtype.casefold():
            return "DUP"
        if "tra" in svtype.casefold():
            return "TRA"
        if "cnv" in svtype.casefold():
            return "CNV"
        return "UNK"
    def print_sv(self):
        print(self.svtype + ": " + self.chrom + " " + str(self.pos) + "(" + str(self.cipos1) +", " + str(self.cipos2) + ")" + " - " + str(self.end) + "(" + str(self.cipos1) +", " + str(self.cipos2) + ")" + " LEN: " + str(self.svlen) + " GT: " + self.gt)
    def checkOverlap(self, sv2):
        if(self.chrom != sv2.chrom):
            return False
        # bear in mind that cipos first coord is negative, hence just addition (example cipos=-10,10)
        minPos1 = self.pos+self.cipos1
        maxPos1 = self.pos+self.cipos2
        minPos2 = sv2.pos+sv2.cipos1
        maxPos2 = sv2.pos+sv2.cipos2

        minEnd1 = self.end+self.ciend1
        maxEnd1 = self.end+self.ciend2
        minEnd2 = sv2.end+self.ciend1
        maxEnd2 = sv2.end+self.ciend2

        max_between = 100
        if(self.svtype == "INS"):
            if(abs(self.pos-sv2.pos) < max_between and abs(self.svlen-sv2.svlen) < max_between):
                return True 
        else:
            if(max(minPos1, minPos2)-max_between <= min(maxPos1, maxPos2)):
                if(max(minEnd1, minEnd2)-max_between <= min(maxEnd1, maxEnd2)):
                    return True
        return False

class SVTool:
    max_conf = 200 # max confidence interval length
    def __init__(self, filename):
        self.tool = filename.split("/")[-1].split(".")[0]
        self.parse_file(filename)
    def parse_file(self, filename):
        self.sv_list = list()
        with open(filename) as file:
            for line in file:
                if not(line.startswith('#')):
                    sv = SVariant(self.tool, line)
                    if(abs(sv.ciend2-sv.ciend1) > self.max_conf or abs(sv.cipos2-sv.cipos1) > self.max_conf):
                        continue
                    #print(self.tool + " | ", end = '')
                    #sv.print_sv()
                    self.sv_list.append(sv)