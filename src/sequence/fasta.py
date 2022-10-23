import re
from Bio import SeqIO


fasta_rules = {
    'default': \
        r'^(?P<id>\S+)(?: (?P<description>.+))?',
    'UniProt': \
        r'^(?P<database>\w+)\|(?P<id>\w+)\|(?P<name>\w+)' + \
        '(?: ' + \
        r'(?P<description>(?:(?![A-Z]{2}=).)+) ' + \
        r'OS=(?P<organism>(?:(?![A-Z]{2}=).)+) ' + \
        r'OX=(?P<organismId>(?:(?![A-Z]{2}=).)+) ' + \
        r'(GN=(?P<gene>(?:(?![A-Z]{2}=).)+) )?' + \
        r'PE=(?P<proteinExistence>(?:(?![A-Z]{2}=).)+) ' + \
        r'SV=(?P<sequenceVersion>(?:(?![A-Z]{2}=).)+)' + \
        ')?'
}


def build_fasta_header_parser(parsing_rule):
    if type(parsing_rule).__name__.find('Pattern') >= 0:
        regex = parsing_rule
    else:
        rule = fasta_rules.get(parsing_rule, None)
        if rule is None:
            raise ValueError('invalid protease: ' + str(parsing_rule))
        regex = re.compile(rule)

    def _parse(header):
        match = regex.search(header)
        if match is not None:
            return match.groupdict()
        else:
            None

    return _parse



def read_fasta(file, parsing_rule='default', return_iterator=False):
    header_parser = build_fasta_header_parser(parsing_rule)

    def _parse(record):
        result = header_parser(record.description)

        if result is None:
            result = {
                'id': record.id,
                'description': record.description
            }

        result['sequence'] = str(record.seq)

        return result

    fasta = SeqIO.parse(file, 'fasta')
    proteins = (_parse(record) for record in fasta)

    if not return_iterator:
        proteins = list(proteins)

    return proteins


