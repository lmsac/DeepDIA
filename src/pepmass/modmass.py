from .pepmass import PeptideMassCalculator
from .modinfo import ModInfo, find_modification, \
    default_fixed_modifications, default_variable_modifications


class ModSite:
    def __init__(self, name, position, site=None):
        self.name = name
        self.position = position
        self.site = site

    def __str__(self):
        s = self.name
        if self.site == 'N-term' or self.site == 'C-term':
            s += '(' + self.site + ')'
        elif self.site is not None:
            s += '(' + self.site + str(self.position) + ')'
        else:
            s += '(' + str(self.position) + ')'
        return s

    @staticmethod
    def n_term(name):
        return ModSite(name, None, 'N-term')

    @staticmethod
    def c_term(name):
        return ModSite(name, None, 'C-term')

    @staticmethod
    def from_dict(d):
        return ModSite(
            d['name'],
            d.get('position', None),
            d.get('site', None)
        )
    

class ModifiedPeptideMassCalculator(PeptideMassCalculator):
    def __init__(self,
                 fixed_modifications=None,
                 variable_modifications=None,
                 **kwargs):
        super(ModifiedPeptideMassCalculator, self) \
            .__init__(**kwargs)

        if fixed_modifications is None:
            fixed_modifications = default_fixed_modifications()
        
        self.fixed_modifications = fixed_modifications

        if variable_modifications is None:
            variable_modifications = default_variable_modifications()
        
        self.variable_modifications = variable_modifications


    def find_var_mod(self, sequence, modification):
        if isinstance(modification, dict):
            modification = ModSite.from_dict(modification)

        site = sequence[modification.position - 1] \
            if isinstance(modification.position, int) else None

        if isinstance(modification.site, str):
            if site is None:
                site = modification.site
            elif modification.site != site:
                raise ValueError(
                    'sequence and modification site not match: ' + \
                    sequence + ', ' + str(modification)
                )

        mod = find_modification(
            self.variable_modifications,
            modification.name, site
        )
        if mod is not None:
            return mod

        mod = find_modification(
            self.fixed_modifications,
            modification.name, site
        )
        if mod is not None:
            # is a fixed modification
            return None
        else:
            raise ValueError(
                'unknown modification: ' + \
                modification.name + '[' + site + ']'
            )


    def mw(self, sequence, modification=None, **kwargs):
        def _fixed_mod_count(sequence, mod, site=None):
            if isinstance(mod, ModInfo):
                site = mod.site
            if site == 'N-term' or site == 'C-term':
                return 1
            elif isinstance(site, str):
                return sequence.count(site)
            elif isinstance(site, list) or isinstance(site, tuple):
                return sum((
                    _fixed_mod_count(sequence, None, site=site)
                    for s in site
                ))
            else:
                raise TypeError('invalid site: ' + str(type(site)))

        def _fixed_mod_mass(sequence):
            return sum((
                mod.delta_mass * _fixed_mod_count(sequence, mod)
                for mod in self.fixed_modifications
            ))

        def _var_mod_mass(sequence, modification):
            if modification is None:
                return 0.0

            if isinstance(modification, dict):
                modification = ModSite.from_dict(modification)

            if isinstance(modification, ModSite):
                mod = self.find_var_mod(sequence, modification)
                if mod is None:
                    return 0.0
                else:
                    return mod.delta_mass
            elif isinstance(modification, list):
                return sum((
                    _var_mod_mass(sequence, mod)
                    for mod in modification
                ))
            else:
                raise TypeError('invalid modification: ' + \
                    str(type(modification)))

        result = super(ModifiedPeptideMassCalculator, self) \
            .mw(sequence=sequence, **kwargs)

        result += _fixed_mod_mass(sequence)
        result += _var_mod_mass(sequence, modification)
        return result


    def fragment_mod_count(self, sequence, modification=None,
                           fragment_type = 'b'):
        def _fixed_mod_count(sequence):
            def __fixed_mod_count(sequence, mod):
                if mod.site == 'N-term':
                    return [1] * len(sequence)
                elif mod.site == 'C-term':
                    return [0] * (len(sequence) - 1) + [1]

                cumsum_count = []
                sum_count = 0
                for aa in sequence:
                    if aa == mod.site or aa in mod.site:
                        sum_count += 1
                    cumsum_count.append(sum_count)
                return cumsum_count

            return {
                x[0]: x[1]
                for x in (
                    (mod, __fixed_mod_count(sequence, mod))
                    for mod in self.fixed_modifications
                )
                if sum(x[1]) > 0
            }

        def _var_mod_count(sequence, modification):
            count_dict = dict()
            if isinstance(modification, dict):
                modification = ModSite.from_dict(modification)
            if isinstance(modification, ModSite):
                modification = [modification]

            if isinstance(modification, list):
                for m in modification:
                    if isinstance(m, dict):
                        m = ModSite.from_dict(m)
                    elif not isinstance(m, ModSite):
                        raise TypeError('invalid modification: ' + \
                            str(type(m)))
                    mod = self.find_var_mod(sequence, m)
                    if mod is None:
                        continue
                    cumsum_count = count_dict.get(mod, None)
                    if cumsum_count is None:
                        cumsum_count = [0] * (len(sequence))
                        count_dict[mod] = cumsum_count
                    if m.site == 'N-term':
                        for i in range(len(cumsum_count)):
                            cumsum_count[i] += 1
                    elif m.site == 'C-term':
                        cumsum_count[-1] += 1
                    else:
                        for i in range(m.position - 1, len(cumsum_count)):
                            cumsum_count[i] += 1
            elif modification is not None:
                raise TypeError('invalid modification: ' + \
                    str(type(modification)))
            return count_dict


        def _fragment_mod_count_type(mod_count, fragment_type,
                                     allow_list=False):
            if isinstance(fragment_type, str):
                frag_type = self.fragment_type(fragment_type)
                if frag_type.n_term:
                    return {
                        k: v[:-1]
                        for k, v in mod_count.items()
                    }
                else:
                    return {
                        k: [v[-1] - x for x in reversed(v[:-1])]
                        for k, v in mod_count.items()
                    }

            elif allow_list and \
                (isinstance(fragment_type, list) or \
                isinstance(fragment_type, tuple)):
                result = [
                    {
                        'fragment_type': t,
                        'mod_count': _fragment_mod_count_type(mod_count, t)
                    }
                    for t in fragment_type
                ]
                return result
            else:
                raise TypeError('invalid fragment_type: ' + \
                    str(type(fragment_type)))

        fixed_mod_count = _fixed_mod_count(sequence)
        var_mod_count = _var_mod_count(sequence, modification)

        mod_count = dict()
        mod_count.update(fixed_mod_count)
        mod_count.update(var_mod_count)

        result = _fragment_mod_count_type(
            mod_count, fragment_type,
            allow_list=True
        )
        return result


    def fragment_neutral_mw(self, sequence, modification=None,
                            fragment_type='b', loss=None,
                            **kwargs):
        def _fragment_mw_mod(fragment_mw, mod_count):
            if (mod_count is None):
                return fragment_mw

            if not isinstance(fragment_mw, list):
                raise TypeError('invalid fragment_mw: ' + \
                    str(type(fragment_mw)))
            if len(fragment_mw) == 0:
                return fragment_mw

            if isinstance(fragment_mw[0], float):
                if not isinstance(mod_count, dict):
                    raise TypeError('invalid mod_count: ' + \
                        str(type(mod_count)))

                result = fragment_mw.copy()
                for k, v in mod_count.items():
                    for i, x in enumerate(v):
                        result[i] += k.delta_mass * x
                return result

            elif isinstance(fragment_mw[0], dict):
                def update(d, **kwargs):
                    d.update(kwargs)
                    return d

                if isinstance(mod_count, dict):
                    return [
                        update(
                            x.copy(),
                            fragment_mw=_fragment_mw_mod(
                                x['fragment_mw'], mod_count
                            )
                        )
                        for x in fragment_mw
                    ]
                elif isinstance(mod_count, list):
                    def match_record(x, lst):
                        for i, y in enumerate(lst):
                            if all(y.get(k, v) == v for k, v in x.items()):
                                return i, y
                        return None, None

                    result = [
                        update(
                            t[0].copy(),
                            fragment_mw=_fragment_mw_mod(
                                t[0]['fragment_mw'],
                                t[1]['mod_count']
                            )
                        )
                        for t in
                        (
                            (x, match_record(x, mod_count)[1])
                            for x in fragment_mw
                        )
                        if t[1] is not None
                    ]
                    return result
                else:
                    raise TypeError('invalid mod_count: ' + \
                        str(type(mod_count)))

            else:
                raise TypeError('invalid fragment_mw: ' + \
                    str(type(fragment_mw)))


        def _fragment_loss_mass(mod_count, loss_count=None):
            if isinstance(mod_count, dict):
                loss_mass_list = []

                for k, v in mod_count.items():
                    if k.loss_mass is None:
                        continue
                    n_loss = None
                    if isinstance(loss_count, dict):
                        if k not in loss_count:
                            continue
                        n_loss = loss_count[k]
                    if n_loss is None:
                        from itertools import count
                        n_loss = count(1)

                    for n in n_loss:
                        if n > max(v):
                            break
                        loss_mass_list.append((
                            (str(n) + '*' if n > 1 else '') + k.loss_name,
                            [k],
                            [
                                n * k.loss_mass if x >= n else None
                                for j, x in enumerate(v)
                            ]
                        ))
                        for i in range(len(loss_mass_list)):
                            if k in loss_mass_list[i][1]:
                                continue
                            loss_mass_list.append((
                                loss_mass_list[i][0] + '+' + \
                                (str(n) + '*' if n > 1 else '') + \
                                k.loss_name,
                                loss_mass_list[i][1] + [k],
                                [
                                    loss_mass_list[i][2][j] + \
                                    n * k.loss_mass \
                                    if x >= n and \
                                    loss_mass_list[i][2][j] is not None \
                                    else None
                                    for j, x in enumerate(v)
                                ]
                            ))

                result = {
                    t[0]: t[2]
                    for t in loss_mass_list
                }
                return result

            elif isinstance(mod_count, list):
                def update(d, **kwargs):
                    d.update(kwargs)
                    return d
                def remove(d, key):
                    d.pop(key)
                    return d

                return [
                    update(
                        remove(x.copy(), 'mod_count'),
                        loss_mass=_fragment_loss_mass(
                            x['mod_count'],
                            loss_count
                        )
                    )
                    for x in mod_count
                ]
            else:
                raise TypeError('invalid mod_count: ' + \
                    str(type(mod_count)))


        def _fragment_mw_mod_loss(fragment_mw, loss_mass):
            if (loss_mass is None):
                return fragment_mw

            if not isinstance(fragment_mw, list):
                raise TypeError('invalid fragment_mw: ' + \
                    str(type(fragment_mw)))
            if len(fragment_mw) == 0:
                return fragment_mw

            if isinstance(fragment_mw[0], float):
                if not isinstance(loss_mass, dict):
                    raise TypeError('invalid loss_mass: ' + \
                        str(type(loss_mass)))
                return [
                    {
                        'loss': k,
                        'fragment_mw': [
                            fragment_mw[i] - x if x is not None else None
                            for i, x in enumerate(v)
                        ]
                    }
                    for k, v in loss_mass.items()
                ]

            elif isinstance(fragment_mw[0], dict):
                def update(d, **kwargs):
                    d.update(kwargs)
                    return d

                if isinstance(loss_mass, dict):
                    result = []
                    for x in fragment_mw:
                        loss = x.get('loss', None)
                        if loss is None or loss == 'noloss':
                            loss = ''
                        else:
                            loss = '+' + loss
                        result.extend(
                            update(
                                x.copy(),
                                fragment_mw=r['fragment_mw'],
                                loss=r['loss'] + loss
                            )
                            for r in _fragment_mw_mod_loss(
                                x['fragment_mw'], loss_mass
                            )
                        )
                    return result
                elif isinstance(loss_mass, list):
                    def match_record(x, lst):
                        for i, y in enumerate(lst):
                            if all(y.get(k, v) == v for k, v in x.items()):
                                return i, y
                        return None, None

                    result = []
                    for x in fragment_mw:
                        loss = x.get('loss', None)
                        if loss is None or loss == 'noloss':
                            loss = ''
                        else:
                            loss = '+' + loss
                        l = match_record(x, loss_mass)[1]
                        if l is not None:
                            result.extend(
                                update(
                                    x.copy(),
                                    fragment_mw=r['fragment_mw'],
                                    loss=r['loss'] + loss
                                )
                                for r in _fragment_mw_mod_loss(
                                    x['fragment_mw'], l['loss_mass']
                                )
                            )
                    return result

                else:
                    raise TypeError('invalid mod_count: ' + \
                        str(type(mod_count)))
            else:
                raise TypeError('invalid fragment_mw: ' + \
                    str(type(fragment_mw)))


        def _parse_loss(loss):
            def __parse_loss_single(loss):
                if loss == 'modloss':
                    return None, 'any'
                if loss in self.neutral_losses:
                    return loss, None
                for m in self.fixed_modifications + self.variable_modifications:
                    if m.loss_name == loss:
                        return None, loss
                return None, None

            def __parse_loss_term(loss):
                c, m = __parse_loss_single(loss)
                if c is not None or m is not None:
                    return c, m
                t = loss.split('*', 1)
                if len(t) == 2:
                    if t[0].isdigit():
                        n = int(t[0])
                        if n > 0:
                            c, m = __parse_loss_single(t[1])
                            if m is not None and m != 'any':
                                return None, {m: n}
                return None, None

            def __parse_loss_sum(loss):
                c, m = __parse_loss_term(loss)
                if c is not None or m is not None:
                    return c, m
                common_loss = None
                mod_loss = None
                for s in loss.split('+'):
                    c, m = __parse_loss_term(s)
                    if c is not None:
                        if common_loss != None:
                            return None, None
                        common_loss = c
                    elif m is not None:
                        if mod_loss == 'any':
                            return None, None
                        if m == 'any':
                            if mod_loss is not None:
                                return None, None
                            mod_loss = 'any'
                            continue
                        if mod_loss is None:
                            mod_loss = m
                            continue
                        if isinstance(mod_loss, str):
                            mod_loss = {mod_loss: 1}
                        if isinstance(m, str):
                            mod_loss[m] = mod_loss.get(m, 0) + 1
                        elif isinstance(m, dict):
                            for k, v in m:
                                mod_loss[k] = mod_loss.get(k, 0) + v
                        else:
                            return None, None
                    else:
                        return None, None
                return common_loss, mod_loss

            if loss is None:
                return 'noloss', None
            if loss == '' or loss == 'noloss' or loss == 'None':
                return 'noloss', None

            if isinstance(loss, str):
                c, m = __parse_loss_sum(loss)
                if c is not None or m is not None:
                    return c, m
                else:
                    raise ValueError('invalid loss: ' + loss)

            if isinstance(loss, list) or \
                isinstance(loss, tuple):
                return [
                    _parse_loss(l)
                    for l in loss
                ]
            else:
                raise TypeError('invalid loss: ', str(type(loss)))


        def _filter_loss_mass(loss_mass, loss):
            def __get_loss_id(loss):
                if isinstance(loss, tuple):
                    mod_loss = loss[1]
                    if mod_loss == 'any':
                        return mod_loss
                    if isinstance(mod_loss, dict):
                        return '+'.join(
                            (str(v) + '*' if v > 1 else '') + k
                            for k, v in mod_loss.items()
                        )
                    else:
                        return mod_loss
                elif isinstance(loss, list):
                    result = set()
                    for l in loss:
                        r = __get_loss_id(l)
                        if r == 'any':
                            return 'any'
                        if r is not None:
                            result.add(r)
                    return result
                else:
                    return loss

            loss = __get_loss_id(loss)

            if loss == 'any':
                return loss_mass
            if loss is None:
                return None

            if isinstance(loss_mass, dict):
                if isinstance(loss, str):
                    value = loss_mass.get(loss, None)
                    if value is not None:
                        return {
                            loss: value
                        }
                    else:
                        return {}
                else:
                    return {
                        k: v for k, v in loss_mass.items()
                        if k in loss
                    }
            elif isinstance(loss_mass, list):
                def update(d, **kwargs):
                    d.update(kwargs)
                    return d
                return [
                    update(
                        x.copy(),
                        loss_mass=_filter_loss_mass(x['loss_mass'], loss)
                    )
                    for x in loss_mass
                ]
            else:
                raise TypeError('invalid loss_mass: ', str(type(loss_mass)))


        def _filter_fragment_mw_loss(fragment_mw, loss):
            def __get_loss_id(loss):
                if isinstance(loss, tuple):
                    common_loss = loss[0]
                    mod_loss = loss[1]
                    if mod_loss is None:
                        return common_loss
                    if mod_loss == 'any':
                        return None
                    if isinstance(mod_loss, dict):
                        mod_loss = '+'.join(
                            (str(v) + '*' if v > 1 else '') + k
                            for k, v in mod_loss.items()
                        )
                    return mod_loss + \
                        ('+' + common_loss if common_loss is not None else '')
                elif isinstance(loss, list):
                    result = set()
                    for l in loss:
                        r = __get_loss_id(l)
                        if r is not None:
                            result.add(r)
                    return result
                else:
                    return loss

            def __get_common_loss_id_any_mod_loss(loss):
                if isinstance(loss, tuple):
                    mod_loss = loss[1]
                    if mod_loss == 'any':
                        common_loss = loss[0]
                        if common_loss is None:
                            common_loss = 'noloss'
                        return common_loss
                    return None
                elif isinstance(loss, list):
                    result = set()
                    for l in loss:
                        r = __get_common_loss_id_any_mod_loss(l)
                        if r is not None:
                            result.add(r)
                    return result
                else:
                    return loss

            common_loss_any_mod_loss = __get_common_loss_id_any_mod_loss(loss)
            loss = __get_loss_id(loss)

            if isinstance(fragment_mw, list) and \
                len(fragment_mw) > 0 and \
                isinstance(fragment_mw[0], dict):
                if isinstance(loss, str):
                    loss_func = lambda x: x == loss
                elif loss is not None:
                    loss_func = lambda x: x in loss
                else:
                    loss_func = lambda x: False
                if isinstance(common_loss_any_mod_loss, str):
                    common_loss_any_func = lambda x: \
                        x.endswith('+' + common_loss_any_mod_loss)
                elif common_loss_any_mod_loss is not None:
                    if 'noloss' in common_loss_any_mod_loss:
                        common_loss_any_func = lambda x: \
                            not any(map(lambda l:
                                x.get('loss', '').endswith('+' + l),
                                self.neutral_losses.keys()))
                    else:
                        common_loss_any_func = lambda x: \
                            any(map(
                                lambda l: x.endswith('+' + l),
                                common_loss_any_mod_loss)
                            )
                else:
                    common_loss_any_func = lambda x: False

                return [
                    x for x in fragment_mw
                    if loss_func(x.get('loss', '')) or \
                        common_loss_any_func(x.get('loss', ''))
                ]
            else:
                return fragment_mw

        loss = _parse_loss(loss)
        if isinstance(loss, list):
            has_mod_loss = any(map(lambda t: t[1] is not None, loss))
        else:
            has_mod_loss = loss[1] is not None
        if isinstance(loss, list):
            common_loss = list(set(map(lambda t: t[0], loss)))
        elif has_mod_loss:
            common_loss = [loss[0]]
        else:
            common_loss = loss[0]


        fragment_mw = super(ModifiedPeptideMassCalculator, self) \
            .fragment_neutral_mw(
                sequence=sequence,
                fragment_type=fragment_type,
                loss=common_loss,
                **kwargs
            )

        mod_count = self.fragment_mod_count(
            sequence=sequence, 
            modification=modification,
            fragment_type=fragment_type
        )
        fragment_mw = _fragment_mw_mod(fragment_mw, mod_count)
        
        if has_mod_loss:
            loss_mass = _fragment_loss_mass(mod_count, None)
            loss_mass = _filter_loss_mass(loss_mass, loss)
            if loss_mass is not None:
                fragment_mw.extend(
                    _fragment_mw_mod_loss(fragment_mw, loss_mass)
                )
                fragment_mw = _filter_fragment_mw_loss(fragment_mw, loss)

        return fragment_mw



if __name__ == '__main__':
    pep_calc = ModifiedPeptideMassCalculator();

    print(pep_calc.fragment_neutral_mw(
        'LCISWYDNEFGTYSNR',
        modification=[
            {'name': 'Phospho', 'position': 4},
            {'name': 'Acetyl', 'site':'N-term'},
            {'name': 'Carbamidomethyl', 'position': 2},
            {'name': 'Phospho', 'position': 6},
            {'name': 'Phospho', 'position': 12}
        ],
        fragment_type='b',
        #loss = 'H3PO4+H2O',
        loss = ['NH3', 'H3PO4+H2O', 'modloss+NH3', 'noloss']
    ))

