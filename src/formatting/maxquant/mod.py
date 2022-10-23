import re

from pepmass.modinfo import find_modification, \
    default_fixed_modifications, default_variable_modifications


MaxQuant_modification_abbreviation = {
    'ox': 'Oxidation', 
    'ac': 'Acetyl'
}


def parse_modification_MaxQuant(
    modified_sequence, 
    mod_abbr=MaxQuant_modification_abbreviation):
    mods = []

    matches = re.findall(
        '([_A-Z])(?:\\(([^\\(\\)]+)(?:\\)))?',
        modified_sequence
    )

    if matches[0][1] == '_':
        if matches[0][2]:
            mods.append({
                'name': mod_abbr.get(matches[0][1], matches[0][1]),
                'site': 'N-term',
                'position': None
            })
        matches.pop(0)

    for i, m in enumerate(matches):
        if m[1]:
            if m[0] != '_':
                mods.append({
                    'name': mod_abbr.get(m[1], m[1]),
                    'site': m[0],
                    'position': i
                })
            elif i == len(matches) - 1:
                mods.append({
                    'name': mod_abbr.get(m[1], m[1]),
                    'site': 'C-term',
                    'position': None
                })

    return mods or None


def stringify_modinfo_MaxQuant(modinfo):
    s = modinfo.name
    if isinstance(modinfo.site, str):
        s += ' (' + modinfo.site + ')'
    elif isinstance(modinfo.site, list):
        s += ' (' + ''.join(modinfo.site) + ')'
    return s


def to_modified_sequence_MaxQuant(
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
                    s = stringify_modinfo_MaxQuant(modinfo)

            if s is None and fixed_modifications is not None:
                modinfo = find_modification(
                    fixed_modifications,
                    name=mod['name'],
                    site=mod['site'] or aa[mod['position']]
                )
                if modinfo is not None:
                    continue

            if s is None:
                s = mod['name'] + '(' + (mod['site'] or aa[mod['position']]) + ')'

            if mod.get('site', None) == 'N-term':
                aa[0] += '('+ s + ')'
            elif mod.get('site', None) == 'C-term':
                aa[-1] += '('+ s + ')'
            else:
                aa[mod['position']] += '('+ s + ')'

    return ''.join(aa)

