import argparse, sys
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
import cyvcf2

class UpdateInfo(object):
    def __init__(self, vcf_filename, output_filename):
        self.vcf_filename = vcf_filename
        self.output_filename = output_filename

    def calc_msq(self, var):
            # Below is what was in vcfpaste, but what if multiple ALTs?
            # Do we only expect 1 ALT per line?
            # NOTE SQ is defined as: 'Phred-scaled probability that this site is variant (non-reference in this sample'
            # Likely want average sample quality across all non-0/0 genotypes rather than just those containing 1
            gt = [var.genotype(s).get_format('GT') for s in var.sample_list]
            positive_gt = filter(lambda x: x == '0/1' or x == '1/1', gt)
            num_pos = len(positive_gt)
            sum_sq = 0.0
            try:
                sum_sq += sum([float(var.genotype(s).get_format('SQ')) for s in var.sample_list \
                if var.genotype(s).get_format('GT') == '1/1' or var.genotype(s).get_format('GT') == '0/1'])
            except ValueError:
                sum_sq += 0
            if num_pos > 0:
                msq = '%0.2f' % (sum_sq / num_pos)
            else:
                msq = '.'
            return msq

    def execute(self):
        vcf = cyvcf2.VCF(self.vcf_filename)
        vcf.add_info_to_header({
            'ID' : 'AF', 
            'Number' : 'A', 
            'Type' : 'Float', 
            'Description' : 'Allele Frequency, for each ALT allele, in the same order as listed'
            })
        vcf.add_info_to_header({
            'ID' : 'NSAMP',
            'Number' : '1', 
            'Type' : 'Integer',
            'Description' : 'Number of samples with non-reference genotypes'
            })
        vcf.add_info_to_header({
            'ID' : 'MSQ', 
            'Number' : '1', 
            'Type' : 'Float', 
            'Description' : 'Mean sample quality of positively genotyped samples'
            })

        output = cyvcf2.Writer(self.output_filename, vcf)

        for variant in vcf:
            af = variant.aaf
            nsamp = variant.num_het + variant.num_hom_alt
            variant.INFO['AF'] = '%0.4f' % af
            variant.INFO['NSAMP'] = str(nsamp)
            msq = None
            try:
                msq = '%0.2f' % variant.msq
            except ZeroDivisionError:
                msq = '.'
            variant.INFO['MSQ'] = msq
            output.write_record(variant)

def description():
    return 'Add allele frequency information to a VCF file'

def add_arguments_to_parser(parser):
    parser.add_argument(metavar='vcf', dest='input_vcf', nargs='?', default=None, help='VCF input')
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description=description())
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    handle = None
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input_vcf = '-'
    updater = UpdateInfo(args.input_vcf, '-')
    updater.execute()

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))
