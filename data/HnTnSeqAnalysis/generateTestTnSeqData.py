from optparse import OptionParser


def main(options, args):



if __name__ == '__main__':
   parser = OptionParser()
   parser.add_option("-o", "--outDir", dest="filename", default='./testTnSeqData', help="directory for output", metavar="FILE")
   parser.add_option("-o", "--outDir", dest="filename", default='./testTnSeqData', help="directory for output", metavar="FILE")

   (options, args) = parser.parse_args()
   main(options, args)