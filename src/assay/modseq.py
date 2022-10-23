import re


def stringify_modification(modification):
    if modification is None or modification == 'null' or \
        modification == 'None':
        return None

    def _to_str(m):
        r = m['name'] + '('
        if m['position'] is not None:
            r += str(m['position'])
        if m['site'] is not None:
            r += str(m['site'])
        r += ')'
        return r

    return ';'.join ((
        _to_str(m)
        for m in modification
    ))


def parse_modification(s):
    if s is None or \
        s == 'null' or \
        s == 'None':
        return None

    def _from_str(s):
        mat = re.match('^(.+)\\(([0-9]+)?([A-Z]|N-term|C-term)?\\)$', s)
        name = mat.group(1)
        pos = mat.group(2)
        if pos is not None:
            pos = int(pos)
        site = mat.group(3)
        return {
            'name': name,
            'position': pos,
            'site': site
        }

    return [
        _from_str(x)
        for x in s.split(';')
        if x != ''
    ]


def check_modifications(modifications,
                        action='if_any',
                        name_list=None):
        if action == 'if_any' or \
            action == 'include':
            if modifications is None or len(modifications) == 0:
                return False
            elif name_list is None or \
                len(name_list) == 0:
                return True
            elif not any(
                m['name'] in name_list
                for m in modifications
            ):
                return False

        elif action == 'if_exclusive':
            if modifications is None or len(modifications) == 0:
                return False
            elif name_list is None or \
                len(name_list) == 0:
                return True
            elif not any(
                m['name'] not in name_list
                for m in modifications
            ):
                return False

        elif action == 'not_if_any' or \
            action == 'exclude':
            if modifications is None or len(modifications) == 0:
                return True
            elif name_list is None or \
                len(name_list) == 0:
                return False
            elif any(
                m['name'] in name_list
                for m in modifications
            ):
                return False

        elif action == 'not_if_exclusive':
            if modifications is None or len(modifications) == 0:
                return True
            elif name_list is None or \
                len(name_list) == 0:
                return False
            elif any(
                m['name'] not in name_list
                for m in modifications
            ):
                return False

        else:
            raise ValueError('invalid modification action: ' + \
                             str(action))

        return True

