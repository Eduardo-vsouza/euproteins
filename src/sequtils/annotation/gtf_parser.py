import os


class GTFParser(object):
    def __init__(self, gtf):
        self.gtf = gtf

    def parse(self):
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split("\t")
                attrs = cols[8].split(";")
                start = cols[3]
                end = cols[4]
                feature = [2]

