import re

from pepmass.modinfo import find_modification, \
    default_fixed_modifications, default_variable_modifications


def parse_modification_Spectronaut(modified_sequence):
    mods = []

    matches = re.findall(
        '(([_A-Z])(?:\\[([^\\(\\)]+)(?: \\([^\\(\\)]+\\))?\\])?)',
        modified_sequence
    )

    if matches[0][1] == '_':
        if matches[0][2]:
            mods.append({
                'name': matches[0][2],
                'site': 'N-term',
                'position': None
            })
        matches.pop(0)

    for i, m in enumerate(matches):
        if m[2]:
            if m[1] != '_':
                mods.append({
                    'name': m[2],
                    'site': m[1],
                    'position': i
                })
            elif i == len(matches) - 1:
                mods.append({
                    'name': m[2],
                    'site': 'C-term',
                    'position': None
                })

    return mods or None


def stringify_modinfo_Spectronaut(modinfo):
    s = modinfo.name
    if isinstance(modinfo.site, str):
        s += ' (' + modinfo.site + ')'
    elif isinstance(modinfo.site, list):
        s += ' (' + ''.join(modinfo.site) + ')'
    return s


def to_modified_sequence_Spectronaut(
    sequence, modification,
    fixed_modifications=default_fixed_modifications(),
    variable_modifications=default_variable_modifications()):

    if modification is None and fixed_modifications is None:
        return '_' + sequence + '_'

    aa = ['_'] + list(sequence) + ['_']

    if modification is not None:
        for mod in modification:
            if mod is None:
                continue

            s = None
            if variable_modifications is not None:
                modinfo = find_modification(
                    variable_modifications,
                    name=mod['name'],
                    site=mod['site'] or aa[mod['position']]
                )
                if modinfo is not None:
                    s = stringify_modinfo_Spectronaut(modinfo)

            if s is None and fixed_modifications is not None:
                modinfo = find_modification(
                    fixed_modifications,
                    name=mod['name'],
                    site=mod['site'] or aa[mod['position']]
                )
                if modinfo is not None:
                    continue

            if s is None:
                s = mod['name'] + (mod['site'] or aa[mod['position']])

            if mod.get('site', None) == 'N-term':
                aa[0] += '['+ s + ']'
            elif mod.get('site', None) == 'C-term':
                aa[-1] = '['+ s + ']'
            else:
                aa[mod['position']] += '['+ s + ']'

    if fixed_modifications is not None:
        for modinfo in fixed_modifications:
            s = stringify_modinfo_Spectronaut(modinfo)
            if modinfo.site == 'N-term':
                aa[0] += '['+ s + ']'
            elif modinfo.site == 'C-term':
                aa[-1] = '['+ s + ']'
            else:
                for i, x in enumerate(aa):
                    if modinfo.site == x or \
                        (isinstance(modinfo.site, list) and x in modinfo.site):
                        aa[i] += '[' + s + ']'

    return ''.join(aa)

