class FragmentTypeInfo:
    def __init__(self, name, n_term=True, atoms=None):
        self.name = name
        self.n_term = n_term
        self.atoms = atoms

    @staticmethod
    def b():
        return FragmentTypeInfo('b', n_term=True, atoms=None)

    @staticmethod
    def y():
        return FragmentTypeInfo('y', n_term=False, atoms={'H': 2, 'O': 1})

    @staticmethod
    def a():
        return FragmentTypeInfo('a', n_term=True, atoms={'C': -1, 'O': -1})

    @staticmethod
    def x():
        return FragmentTypeInfo('x', n_term=False, atoms={'C': 1, 'O': 2})

    @staticmethod
    def c():
        return FragmentTypeInfo('c', n_term=True, atoms={'N': 1, 'H': 3})

    @staticmethod
    def z():
        return FragmentTypeInfo(
            'z', n_term=False,
            atoms={'H': -1, 'N': -1, 'O': 1}
        )


class PeptideMassCalculator:
    def __init__(self, aa_residues=None, elements=None,
                 fragments=None, neutral_losses=None,
                 **kwargs):
        if aa_residues is None:
            aa_residues = {
                'A': 71.037114,
                'R': 156.101111,
                'N': 114.042927,
                'D': 115.026943,
                'C': 103.009185,
                'E': 129.042593,
                'Q': 128.058578,
                'G': 57.021464,
                'H': 137.058912,
                'I': 113.084064,
                'L': 113.084064,
                'K': 128.094963,
                'M': 131.040485,
                'F': 147.068414,
                'P': 97.052764,
                'S': 87.032028,
                'T': 101.047679,
                'U': 150.95363,
                'W': 186.079313,
                'Y': 163.06332,
                'V': 99.068414
            }
        
        self.aa_residues = aa_residues

        if elements is None:
            elements = {
                'H': 1.007825035,
                'C': 12.0,
                'N': 14.003074,
                'O': 15.99491463,
                'proton': 1.007825035 - 0.000549
            }
            
        self.elements = elements

        if fragments is None:
            fragments = {
                'b': FragmentTypeInfo.b(),
                'y': FragmentTypeInfo.y(),
                'a': FragmentTypeInfo.a(),
                'x': FragmentTypeInfo.x(),
                'c': FragmentTypeInfo.c(),
                'z': FragmentTypeInfo.z()
            }
            
        self.fragments = fragments

        if neutral_losses is None:
            neutral_losses = {
                'NH3': self.element_mass('N') + self.element_mass('H') * 3,
                'H2O': self.element_mass('H') * 2 + self.element_mass('O')
            }
            
        self.neutral_losses = neutral_losses


    def aa_residue_mass(self, aa):
        result = self.aa_residues.get(aa, None)
        if result is None:
            raise ValueError(
                'unknown aa residue: ' + str(aa)
            )
        return result

    def element_mass(self, element):
        result = self.elements.get(element, None)
        if result is None:
            raise ValueError(
                'unknown element: ' + str(element)
            )
        return result

    def fragment_type(self, fragment_type):
        result = self.fragments.get(fragment_type, None)
        if result is None:
            raise ValueError(
                'unknown fragment type: ' + str(fragment_type)
            )
        return result

    def neutral_loss_mass(self, loss):
        result = self.neutral_losses.get(loss, None)
        if result is None:
            raise ValueError(
                'unknown neutral loss: ' + str(loss)
            )
        return result


    def mw(self, sequence, **kwargs):
        if not isinstance(sequence, str):
            raise TypeError('invalid sequence: ' + str(type(sequence)))
        elif len(sequence) == 0:
            raise ValueError('empty sequence')

        aa_residues_mass = sum((self.aa_residue_mass(aa) for aa in sequence))

        mw = aa_residues_mass + \
            self.element_mass('H') * 2 + self.element_mass('O')
        return mw


    def precursor_mz(self, sequence, charge=1, **kwargs):
        def _precursor_mz_from_mw(mw, charge, allow_list=False):
            if isinstance(charge, int):
                if charge == 0:
                    raise ValueError('invalid charge: ' + str(charge))
                else:
                    return (mw + charge * self.element_mass('H')) / \
                        abs(charge)
            elif allow_list and \
                (isinstance(charge, list) or \
                isinstance(charge, tuple) or \
                isinstance(charge, range)):
                result = [
                    {
                        'charge': ch,
                        'precursor_mz': _precursor_mz_from_mw(
                            mw, ch
                        )
                    }
                    for ch in charge
                ]
                return result
            else:
                raise TypeError('invalid charge: ' + str(type(charge)))

        mw = self.mw(sequence=sequence, **kwargs)
        precursor_mz = _precursor_mz_from_mw(mw, charge, allow_list=True)
        return precursor_mz


    def fragment_neutral_mw(self, sequence, fragment_type='b', loss='noloss',
                            **kwargs):
        def _residue_sequence_mass(sequence):
            if not isinstance(sequence, str):
                raise TypeError('invalid sequence: ' + str(type(sequence)))
            elif len(sequence) <= 1:
                raise ValueError('sequence length < 2: ' + sequence)

            cumsum_mw = []
            sum_mw = 0.0
            for aa in sequence:
                sum_mw += self.aa_residue_mass(aa)
                cumsum_mw.append(sum_mw)
            return cumsum_mw

        def _fragment_mw_type(residue_sequence_mass, fragment_type,
                              allow_list=False):
            if not isinstance(residue_sequence_mass, list):
                raise TypeError('invalid residue_sequence_mass: ' + \
                    str(type(residue_sequence_mass)))

            if isinstance(fragment_type, str):
                frag_type = self.fragment_type(fragment_type)
                if frag_type.atoms is None:
                    atom_mass = 0.0
                elif isinstance(frag_type.atoms, dict):
                    atom_mass = sum(
                        self.element_mass(k) * v
                        for k, v in frag_type.atoms.items()
                    )
                else:
                    raise TypeError('invalid frag_type.atoms: ' + \
                        str(type(frag_type.atoms)))
                if frag_type.n_term:
                    return [
                        x + atom_mass
                        for x in residue_sequence_mass[:-1]
                    ]
                else:
                    return [
                        residue_sequence_mass[-1] - x + atom_mass
                        for x in reversed(residue_sequence_mass[:-1])
                    ]

            elif allow_list and \
                (isinstance(fragment_type, list) or \
                isinstance(fragment_type, tuple)):
                result = [
                    {
                        'fragment_type': t,
                        'fragment_mw': _fragment_mw_type(
                            residue_sequence_mass, t
                        )
                    }
                    for t in fragment_type
                ]
                return result
            else:
                raise TypeError('invalid fragment type: ' + \
                    str(type(fragment_type)))


        def _fragment_mw_loss(fragment_mw, loss, allow_list=False):
            if not isinstance(fragment_mw, list):
                raise TypeError('invalid fragment_mw: ' + \
                    str(type(fragment_mw)))
            if len(fragment_mw) == 0:
                return fragment_mw

            if loss is None or loss == '' or \
                loss == 'noloss' or loss == 'None':
                if isinstance(fragment_mw[0], dict):
                    for x in fragment_mw:
                        x['loss'] = 'noloss'
                    return fragment_mw
                else:
                    return fragment_mw
            elif isinstance(loss, str):
                loss_mass = self.neutral_loss_mass(loss)
                if isinstance(fragment_mw[0], dict):
                    def update(d, **kwargs):
                        d.update(kwargs)
                        return d

                    return [
                        update(
                            x.copy(),
                            fragment_mw=_fragment_mw_loss(
                                x['fragment_mw'], loss
                            ),
                            loss=loss
                        )
                        for x in fragment_mw
                    ]
                elif isinstance(fragment_mw[0], float):
                    return [
                        x - loss_mass
                        for x in fragment_mw
                    ]
                else:
                    raise TypeError('invalid fragment_mw: ' + \
                        str(type(fragment_mw)))

            elif allow_list and \
                (isinstance(loss, list) or \
                isinstance(loss, tuple)):
                result = []
                for l in loss:
                    r = _fragment_mw_loss(fragment_mw, l)
                    if len(r) == 0:
                        continue
                    if isinstance(r[0], dict):
                        result.extend(r)
                    else:
                        result.append({
                            'loss': l,
                            'fragment_mw': r
                        })
                return result
            else:
                raise TypeError('invalid loss: ' + str(type(loss)))

        residue_sequence_mass = _residue_sequence_mass(sequence)
        fragment_mw = _fragment_mw_type(
            residue_sequence_mass, fragment_type,
            allow_list=True
        )
        fragment_mw = _fragment_mw_loss(
            fragment_mw, loss,
            allow_list=True
        )
        return fragment_mw


    def fragment_mz(self, sequence, fragment_type='b', charge=1, loss='noloss',
                    **kwargs):
        def _fragment_mz_from_mw(fragment_mw, charge, allow_list=False):
            if not isinstance(fragment_mw, list):
                raise TypeError('invalid fragment_mw: ' + \
                    str(type(fragment_mw)))

            if isinstance(charge, int):
                if charge <= 0:
                    raise ValueError('invalid charge: ' + str(charge))

                if isinstance(fragment_mw[0], dict):
                    def update(d, **kwargs):
                        d.update(kwargs)
                        return d
                    def remove(d, key):
                        d.pop(key)
                        return d

                    return [
                        update(
                            remove(x.copy(), 'fragment_mw'),
                            fragment_mz=_fragment_mz_from_mw(
                                x['fragment_mw'], charge
                            ),
                            charge=charge
                        )
                        for x in fragment_mw
                    ]
                elif isinstance(fragment_mw[0], float) or \
                    fragment_mw[0] is None:
                    return [
                        (x + charge * self.element_mass('proton')) / charge
                        if x is not None else None
                        for x in fragment_mw
                    ]
                else:
                    raise TypeError('invalid fragment_mw: ' + \
                        str(type(fragment_mw)))

            elif allow_list and \
                (isinstance(charge, list) or \
                isinstance(charge, tuple) or \
                isinstance(charge, range)):
                result = []
                for ch in charge:
                    r = _fragment_mz_from_mw(fragment_mw, ch)
                    if len(r) == 0:
                        continue
                    if isinstance(r[0], dict):
                        result.extend(r)
                    else:
                        result.append({
                            'charge': ch,
                            'fragment_mz': r
                        })
                return result
            else:
                raise TypeError('invalid charge: ' + str(type(charge)))

        fragment_mw = self.fragment_neutral_mw(
            sequence=sequence, 
            fragment_type=fragment_type, 
            loss=loss,
            **kwargs
        )
        fragment_mz = _fragment_mz_from_mw(
            fragment_mw, charge,
            allow_list=True
        )
        return fragment_mz



if __name__ == '__main__':
    pep_calc = PeptideMassCalculator();

    print(pep_calc.fragment_mz(
        'PEPTIDE',
        fragment_type=['b', 'y'],
        charge=[1, 2],
        loss = ['noloss', 'NH3', 'H2O']
    ))
