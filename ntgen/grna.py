from ntgen import visualization

TAIL5 = 'AAUU'
STEM1 = 'UCUAC'
STEM2 = 'GUAGA'
LINKER = 'U'
LOOP = 'UCUU'  # As


class gRNA:

    part_names = ['tail5', 'stem1', 'loop', 'stem2', 'linker', 'spacer', 'tail3']

    def __init__(self,
                 tail5: str | None = TAIL5,
                 stem1: str | None = STEM1,
                 loop: str | None = LOOP,
                 stem2: str | None = STEM2,
                 linker: str | None = LINKER,
                 spacer: str | None = None,
                 tail3: str | None = None,
                 stem: str | None = None,
                 structure: str | None = None
                 ):

        if stem is not None:
            stem1, stem2 = stem

        if (stem1 is None and stem2 is not None) or (stem1 is not None and stem2 is None):
            raise ValueError(f'stem1 is {stem1} but stem2 is {stem2}')

        parts = [tail5, stem1, loop, stem2, linker, spacer, tail3]
        self.parts = [p.replace('T', 'U') if p is not None else '' for p in parts]
        self.tail5, self.stem1, self.loop, self.stem2, self.linker, self.spacer, self.tail3 = self.parts

        real_parts = [p for p in parts if p is not None]
        if structure is None:
            self._structure = '.' * len(self.parts)
        elif len(structure) == len(real_parts):
            self._structure = []
            i = 0
            for part in parts:
                if part is not None:
                    self._structure.append(structure[i])
                    i += 1
                else:
                    self.structure.append('-')
        else:
            self._structure = structure

        if len(self._structure) not in [len(self.parts), len(self.sequence)]:
            raise ValueError(f'Length mismatch: structure={len(self._structure)}, parts={len(self.parts)}, sequence={len(self.sequence)}')

    @property
    def sequence(self):
        return ''.join(self.parts)

    @property
    def structure(self):
        if len(self._structure) == len(self.parts):
            return ''.join([c * len(p) for c, p in zip(self._structure, self.parts)])
        else:
            return self._structure

    def get_part_structure(self, part_name):
        idx = 0
        for part_name_, part in zip(self.part_names, self.parts):
            if part_name_ != part_name:
                idx += len(part)
            else:
                return self.structure[idx: idx + len(part)]

    @property
    def colors(self):
        return (['yellowgreen'] * len(self.tail5)
                + ['salmon'] * len(self.stem1)
                + ['goldenrod'] * len(self.loop)
                + ['salmon'] * len(self.stem2)
                + ['darkgray'] * len(self.linker)
                + ['skyblue'] * len(self.spacer)
                + ['yellowgreen'] * len(self.tail3)
                )

    def draw(self, **kwargs):
        visualization.draw(self.structure, sequence=self.sequence, colors=self.colors, **kwargs)
