class Codon:

    codon_dict={'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
                'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
                'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
                'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
                'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
                'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                '---':'-','XXX':'-','???':'?'}
                
    mixture_dict={'W':'AT','R':'AG','K':'GT','Y':'CT',
                  'S':'CG','M':'AC','V':'AGC','H':'ATC',
                  'D':'ATG','B':'TGC','N':'ATGC','-':'-'}
    
    
    def __init__(self, codon_string):
        self.codon = codon_string
    
    @classmethod
    def resolveCodon(cls, codon):
        nonmix = []
        if (codon in cls.codon_dict):
            return[codon]
        elif ((codon.count('-') + codon.count('X')) == 3):
            return ['---']
        elif (1 <= codon.count('-') <= 2) or (1 <= codon.count('X') <= 2):
            return['???']
        for base in codon:
            #Check for mixtures
            if (base in cls.mixture_dict):
                if (not nonmix):
                    nonmix = [x for x in cls.mixture_dict[base]]
                else:
                    nonmix=[x +y for x in nonmix for y in cls.mixture_dict[base]]
            else:
                if (not nonmix):
                    nonmix.append(base)
                else:
                    nonmix = [x + base for x in nonmix]
        return nonmix
