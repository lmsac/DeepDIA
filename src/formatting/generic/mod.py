def stringify_modification(modification):
    if modification is None or modification == 'null':
        return None

    def _to_str(m):
        r = ''
        if m['site'] is not None:
            r += str(m['site'])
        if m['position'] is not None:
            r += str(m['position'])
        r += '(' + m['name'] + ')'
        return r

    return ';'.join ((
        _to_str(m)
        for m in modification
    ))
    

def parse_modification(modification):
    if modification is None or modification == 'null':
        return None

    import re
    def _parse_modification(s):
        mat = re.match('^([A-Z])?([0-9]+)\\((.+)\\)$', s)
        pos = int(mat.group(2))
        name = mat.group(3)
        site = mat.group(1)
        return {
            'name': name,
            'position': pos,
            'site': site
        }

    return [
        _parse_modification(s)
        for s in modification.split(';')
        if s != ''
    ]

