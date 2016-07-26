#!/bin/env python

import sys
import gzip
import argparse
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import svtools.utils as su

def set_from_string(string):
    '''
    Convert a comma separated string to a set of strings
    '''
    return set(string.split(","))

def add_based_on_sname(a, b_set):
    a_set = set_from_string(a.get_info('SNAME'))
    if a_set & b_set: #If they share at least one element in the SNAME field
        return True
    else:
        return False

def load_filter_file(filter_file):
    filter_list = list()

    vcf = Vcf()
    header_lines = list()
    in_header = True
    for line in filter_file:
        if in_header:
            header_lines.append(line)
            if line[0:6] == '#CHROM':
                in_header = False
                vcf.add_header(header_lines)
        else:
            v = line.rstrip().split('\t')
            var = Variant(v, vcf)
            filter_list.append((var.var_id, set_from_string(var.get_info('SNAME'))))
    return filter_list

def sname_filter(input_stream, filter_file, output_stream, complement):
    filter_list = load_filter_file(filter_file)

    vcf = Vcf()
    in_header = True    
    header_lines = list()
    sample_list = None
    for line in input_stream:
        if in_header:
            header_lines.append(line)
            if line[0:6] == '#CHROM':
                in_header = False
                vcf.add_header(header_lines)
                # FIXME This is insufficient for general use
                vcf.add_info('FOUND', '.', 'String', 'Variant id in other file')
                output_stream.write(vcf.get_header() + '\n')
        else:
            v = Variant(line.rstrip().split('\t'), vcf)
            found = False
            for f in filter_list:
                if add_based_on_sname(v, f[1]):
                    found = True
                    value_string = ''
                    try:
                        value_string = ','.join(v.get_info('FOUND'), f[0])
                    except KeyError:
                        value_string = f[0]
                    v.set_info('FOUND', value_string)
                    break
            if found != complement:
                output_stream.write(v.get_var_string() + '\n')

def description():
    return 'look for variants sharing the same original call between two VCF files'

def add_arguments_to_parser(parser):
    parser.add_argument('-v', '--complement', action='store_true', dest='complement', default=False, help='return complement of overlap')
    parser.add_argument("-i", "--input", dest="input", metavar='<VCF>', help="VCF file containing variants to be output")
    parser.add_argument("-f", "--filter", dest="filter_file", metavar='<VCF>', type=argparse.FileType('r'), help="VCF file containing variants used to determine if site should be output")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<VCF>', default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    with su.InputStream(args.input) as stream:
        return sname_filter(stream, args.filter_file, args.output, args.complement)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
