import itertools
import pandas as pd
import numpy as np
import re


protease_rules = {
    'Trypsin': '([KR](?=[^P]))',
    'Trypsin/P': '[KR]',
    'LysC': 'K'
}


def cleave(sequence, protease='Trypsin',
           max_missed_cleavages=0,
           min_peptide_length=7, max_peptide_length=50,
           remove_n_term_methionine=True,
           **kwargs):
    if type(protease).__name__.find('Pattern') >= 0:
        regex = protease
    else:
        rule = protease_rules.get(protease, None)
        if rule is None:
            raise ValueError('invalid protease: ' + str(protease))
        regex = re.compile(rule)

    cleaved_sites = [m.start() + 1 for m in regex.finditer(sequence)]

    if len(cleaved_sites) == 0 or cleaved_sites[-1] != len(sequence):
        cleaved_sites.append(len(sequence))

    if remove_n_term_methionine and sequence[0] == 'M':
        if cleaved_sites[0] != 1:
            cleaved_sites.insert(0, 1)
    else:
        cleaved_sites.insert(0, 0)

    cleaved_sequences = [
        sequence[s:e]
        for s, e in zip(
            cleaved_sites[:-1],
            cleaved_sites[1:]
        )
    ]

    peptides = [
        [
            (cleaved_sites[i], ''.join(cleaved_sequences[i:(i + 1 + n)]))
            for i in range(len(cleaved_sequences) - n)
            if min_peptide_length <= \
                cleaved_sites[i + 1 + n] - cleaved_sites[i] <= \
                max_peptide_length
        ]
        for n in range(max_missed_cleavages + 1)
    ]

    return peptides


def digest(proteins,
           term_window_size=False,
           **kwargs):
    def _build_entry(protein, miss, t):
        r = {
            'protein': protein['id'],
            'sequence': t[1],
            'start': t[0] + 1,
            'end': t[0] + len(t[1]),
            'missedCleavages': miss
        }
        if term_window_size and term_window_size > 0:
            r['nTerminal'] = protein['sequence'] \
                [max(t[0] - term_window_size, 0):t[0]] \
                .rjust(term_window_size, '_')
            r['cTerminal'] = protein['sequence'] \
                [(t[0]+ len(t[1]) + 1): \
                 (t[0]+ len(t[1]) + 1 + term_window_size)] \
                .ljust(term_window_size, '_')
        return r

    def _digest(protein):
        return list(itertools.chain.from_iterable(
            (_build_entry(protein, i, t) for t in l)
            for i, l in enumerate(cleave(protein['sequence'], **kwargs))
        ))

    return pd.DataFrame.from_records(
        itertools.chain.from_iterable(
            _digest(protein)
            for protein in proteins
        )
    )


def filter_peptides_by_mass(peptides,
                            min_peptide_mass=None,
                            max_peptide_mass=None,
                            mass_calculator=None):
    if mass_calculator is None:
        from pepmass import PeptideMassCalculator
        mass_calculator = PeptideMassCalculator()

    peptides = peptides.assign(
        mw=peptides['sequence'].map(lambda s: \
            mass_calculator.mw(s) \
            if all(map(lambda aa: aa in mass_calculator.aa_residues, s)) \
            else np.nan
        )
    )

    if min_peptide_mass is not None:
        peptides = peptides.loc[peptides['mw'] >= min_peptide_mass]
    if max_peptide_mass is not None:
        peptides = peptides.loc[peptides['mw'] <= max_peptide_mass]

    return peptides


def group_duplicated_peptides(peptides):
    agg_dict = {
        k: lambda x: ';'.join(map(str, x)) \
            if len(x) > 1 else x.iloc[0]
        for k in peptides.columns
        if k not in ['sequence', 'missedCleavages']
    }
    if 'missedCleavages' in peptides.columns:
        agg_dict['missedCleavages'] = min
    
    return peptides \
        .groupby(by='sequence', as_index=False) \
        .agg(agg_dict)



if __name__ == '__main__':
    proteins = [
        {
            'sequence': \
            'ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSS'
            'GLYSLSSVVTVPSSNFGTQTYTCNVDHKPSNTKVDKTVERKCCVECPPCPAPPVAGPSVF'
            'LFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTFR'
            'VVSVLTVVHQDWLNGKEYKCKVSNKGLPAPIEKTISKTKGQPREPQVYTLPPSREEMTKN'
            'QVSLTCLVKGFYPSDISVEWESNGQPENNYKTTPPMLDSDGSFFLYSKLTVDKSRWQQGN'
            'VFSCSVMHEALHNHYTQKSLSLSPGK',
            'id': 'P01859'
        }
    ]

    peptides = digest(proteins, protease='Trypsin/P', max_missed_cleavages=1)

    peptides = group_duplicated_peptides(peptides)

    print(peptides)