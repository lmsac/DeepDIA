class ModInfo:
    def __init__(self, name, site, delta_mass, 
                 loss_name=None, 
                 loss_mass=None):
        self.name = name
        self.site = site
        self.delta_mass = delta_mass
        self.loss_name = loss_name
        self.loss_mass = loss_mass

    def __str__(self):
        s = self.name
        if isinstance(self.site, str):
            s += '[' + self.site + ']'
        elif isinstance(self.site, list):
            if all(map(lambda x: len(x) == 1, self.site)):
                s += '[' + ''.join(self.site) + ']'
            else:
                s += '[' + '|'.join(self.site) + ']'
        return s
        

def find_modification(modification_list, name, site):
    mod = next((
        mod for mod in modification_list
        if (mod.name == name) and (site in mod.site)
    ), None)
    return mod


def default_fixed_modifications(): 
    return [
        ModInfo(
            name='Carbamidomethyl',
            site='C',
            delta_mass=57.021464
        )
    ]


def default_variable_modifications():
    return [
        ModInfo(
            name='Acetyl',
            site='N-term',
            delta_mass=42.010565
        ),
        ModInfo(
            name='Oxidation',
            site='M',
            delta_mass=15.994915
        ),
        ModInfo(
            name='Deamidated',
            site='N',
            delta_mass=0.984016
        ),
        ModInfo(
            name='Phospho',
            site=['S', 'T'],
            delta_mass=79.966331,
            loss_name='H3PO4',
            loss_mass=97.976896
        ),
        ModInfo(
            name='Phospho',
            site='Y',
            delta_mass=79.966331
        ),
        ModInfo(
            name='Glu->pyro-Glu',
            site='E',
            delta_mass=-18.010565
        ),
        ModInfo(
            name='Gln->pyro-Glu',
            site='Q',
            delta_mass=-17.026549
        ),
        ModInfo(
            name='Label_18O(2)',
            site='C-term',
            delta_mass=4.008491
        )
    ]

