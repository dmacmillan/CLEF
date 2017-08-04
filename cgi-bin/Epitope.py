class Epitope:

    def __init__(self,
        epitope=None,
        protein=None,
        hlas=None,
        start=None,
        end=None,
        source=None,
        r4=None,
        r2=None
    ):
        self.epitope = epitope
        self.protein = protein
        self.hlas = hlas
        self.start = int(start)
        self.end = int(end)
        self.source = source
        self.r4 = r4
        self.r2 = r2
    
    def __str__(self):
        return ('\t').join([
            (',').join(self.epitope),
            self.protein,
            (',').join(self.hlas),
            str(self.start),
            str(self.end),
            self.source,
            (',').join(self.r4),
            (',').join(self.r2)
        ])
        
    def __repr__(self):
        return 'epitope_object'
        # return ('\t').join([
            # (',').join(self.epitope),
            # self.protein,
            # (',').join(self.hlas),
            # str(self.start),
            # str(self.end),
            # self.source,
            # (',').join(self.r4),
            # (',').join(self.r2)
        # ])

    @staticmethod
    def parseEpitopes(efile, header=True, only_proteins=[], except_proteins=[]):
        results = []
        with open(efile, 'r') as f:
            if header:
                f.readline()
            for line in f:
                line = Epitope(*[x.translate(None,'"') for x in line.strip().split('\t')])
                if line.protein not in only_proteins:
                    continue
                line.epitope = line.epitope.split(',')
                line.hlas = line.hlas.split(',')
                line.r4 = line.r4.split(',')
                line.r2 = line.r2.split(',')
                results.append(line)
        return results