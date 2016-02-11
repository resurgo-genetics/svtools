import re
import time

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        self.reference = ''
        self.sample_list = []
        self.sample_indices = dict()
        self.info_list = []
        self.format_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')

    def add_header(self, header):
        for line in header:
            header_array = line.split('=')
            if header_array[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif header_array[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]
                for i in xrange(0, len(self.sample_list)):
                    self.sample_indices[self.sample_list[i]] = i + 9
            else:
                # Handle if line contains any compatible header lines
                self._add_meta(line, header_array)

    def _add_meta(self, line, header_array):
        a = line[line.find('<')+1:line.rfind('>')]
        r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
        args = [b.split('=')[1] for b in r.findall(a)]
        if header_array[0] == '##INFO':
            self.add_info(*args)
        elif header_array[0] == '##ALT':
            self.add_alt(*args)
        elif header_array[0] == '##FORMAT':
            self.add_format(*args)

    # return the VCF header
    def get_header(self, include_samples=True):
        if include_samples:
            header = '\n'.join(['##fileformat=' + self.file_format,
                                '##fileDate=' + time.strftime('%Y%m%d'),
                                '##reference=' + self.reference] + \
                               [i.hstring for i in self.info_list] + \
                               [a.hstring for a in self.alt_list] + \
                               [f.hstring for f in self.format_list] + \
                               ['\t'.join([
                                   '#CHROM',
                                   'POS',
                                   'ID',
                                   'REF',
                                   'ALT',
                                   'QUAL',
                                   'FILTER',
                                   'INFO',
                                   'FORMAT'] + \
                                          self.sample_list
                                      )])
        else:
            header = '\n'.join(['##fileformat=' + self.file_format,
                                '##fileDate=' + time.strftime('%Y%m%d'),
                                '##reference=' + self.reference] + \
                               [i.hstring for i in self.info_list] + \
                               [a.hstring for a in self.alt_list] + \
                               [f.hstring for f in self.format_list] + \
                               ['\t'.join([
                                   '#CHROM',
                                   'POS',
                                   'ID',
                                   'REF',
                                   'ALT',
                                   'QUAL',
                                   'FILTER',
                                   'INFO']
                                          )])
        return header

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def add_info_after(self, insert_id, id, number, type, desc):
        for i in range(0, len(self.info_list)):
            if insert_id == self.info_list[i].id and (id not in [j.id for j in self.info_list]):
                inf = self.Info(id, number, type, desc)
                self.info_list.insert(i + 1, inf)
                return

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_sample(self, name):
        self.sample_list.append(name)
        # XXX Probably slow. Optimize if we start adding lots of samples ever
        self.sample_indices[name] = self.sample_list.index(name) + 9

    # get the VCF column index of a sample
    # NOTE: this is zero-based, like python arrays
    def sample_to_col(self, sample):
        return self.sample_indices[sample]

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'
        def __eq__(self, other):
            return self.hstring == other.hstring

