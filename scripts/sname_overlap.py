import sys
import gzip
import argparse
from svtools.vcf.file import Vcf
from svtools.bedpe import Bedpe

def set_from_string(string):
    '''
    Convert a comma separated string to a set of strings
    '''
    return set(string.split(","))

def info_tag(bedpe, tag):
    '''
    Return contents of a tag in the INFO field.
    Not general for boolean tags
    Should move to bedpe class or a utility class
    '''
    result = re.split('=', ''.join(filter(lambda x: tag + '=' in x, bedpe.info.split(';'))))[1]
    return result

def add_based_on_sname(a_bedpe, b_bedpe):
    a_set = set_from_string(info_tag(a_bedpe, 'SNAME'))
    b_set = set_from_string(info_tag(b_bedpe, 'SNAME'))
    if a_set & b_set: #If they share at least one element in the SNAME field
        a_bedpe.cohort_vars[b_bedpe.name] = b_bedpe.af
        return True
    else:
        return False

def search(aFile, bFile, bedpe_out, pass_prefix, complement):
    bList = list()
    headerObj = Vcf()
        
    if bFile == "stdin":
        bData = sys.stdin
    elif bFile.endswith('.gz'):
        bData = gzip.open(bFile, 'rb')
    else:
        bData = open(bFile, 'r')
    for bLine in bData:
        if bLine.startswith(pass_prefix):
            continue
        bentry = Bedpe(bLine.rstrip().split('\t'))
        bList.append(bentry)
    
    if aFile == "stdin":
        aData = sys.stdin
    elif aFile.endswith('.gz'):
        aData = gzip.open(aFile, 'rb')
    else:
        aData = open(aFile, 'r')

    in_header = True    
    header_lines = []
    sample_list = None
    for aLine in aData:
        if pass_prefix is not None and aLine.startswith(pass_prefix):
            if aLine[0] == '#' and aLine[1] != '#':
                sample_list = aLine.rstrip().split('\t', 14)[-1]
            else:
                header_lines.append(aLine)
            continue
        else:
            if in_header == True:
                headerObj.add_header(header_lines)

                header = headerObj.get_header()
                bedpe_out.write(header[:header.rfind('\n')] + '\n')                
                if len(sample_list) > 0:
                    bedpe_out.write('\t'.join(['#CHROM_A',
                                               'START_A',
                                               'END_A',
                                               'CHROM_B',
                                               'START_B',
                                               'END_B',
                                               'ID',
                                               'QUAL',
                                               'STRAND_A',
                                               'STRAND_B',
                                               'TYPE',
                                               'FILTER',
                                               'INFO_A','INFO_B',
                                               sample_list]
                                             ) + '\n')
                else:
                    bedpe_out.write('\t'.join(['#CHROM_A',
                                               'START_A',
                                               'END_A',
                                               'CHROM_B',
                                               'START_B',
                                               'END_B',
                                               'ID',
                                               'QUAL',
                                               'STRAND_A',
                                               'STRAND_B',
                                               'TYPE',
                                               'FILTER',
                                               'INFO_A','INFO_B']
                                              ) + '\n')
                in_header=False
            a = Bedpe(aLine.rstrip().split('\t'))
            for b in bList:
                if add_based_on_sname(a, b) != complement:
                    bedpe_out.write(aLine)
                    break

def description():
    return 'look for variants common between two BEDPE files'

def add_arguments_to_parser(parser):
    parser.add_argument('-v', '--complement', action='store_true', dest='complement', default=False, help='return complement of overlap')
    parser.add_argument("-a", "--aFile", dest="aFile", metavar='<BEDPE>', help="pruned, merged BEDPE (A file) or standard input (-a stdin).")
    parser.add_argument("-b", "--bFile", dest="bFile", metavar='<BEDPE>', help="pruned merged BEDPE (B file) (-b stdin). For pruning use svtools prune")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<BEDPE>', default=sys.stdout, help='output BEDPE to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    pass_prefix = "#"
    if args.aFile == None:
        if sys.stdin.isatty():
            sys.stderr.write('Please stream in input to this command or specify the file to read\n')
            sys.exit(1)
        else:
            args.aFile = sys.stdin

    try:
        search(args.aFile, args.bFile, args.output, pass_prefix, complement)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n");

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
