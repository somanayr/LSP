## A class to read FASTA files and act as a container for sequence data.
## @author Sara Thiebaud

class FASTAFile (object):
    """A class to represent sequence data for both proteins and nucleic acids."""

    def __iter__(self):
        self._read()
        return iter(self.records)
    
    def __init__(self, fnp):
        """Create a new instance associated with the specified file.
        You may give either a filename (string) or an open file.
        """
        if isinstance(fnp, str):
            self.fp = open(fnp, 'r')
        elif isinstance(fnp, file):
            self.fp = fnp
        else:
            raise TypeError("Parameter must be a filename or file object")

        self.records = []
        self._read()
    
    def _read(self):
        ln = self.fp.readline()
        while(ln != ''):
            ln = ln.strip()
            if ln[:1] == '>':
                name = ln[1:]
                seq = ''
                ln = self.fp.readline()
                while ln != '' and ln[:1] != '>':
                    seq += ln.strip()
                    ln = self.fp.readline()
                self.records.append(FastaRec(name, seq))
            else:
                ln = self.fp.readline()
    
    def __len__(self):
        return len(self.records)
    
    def __getitem__(self, pos):
        return self.records.__getitem__(pos)
        
class FastaRec:
    def __init__(self, n, s):
        self.name = n
        self.seq = s.upper()
        
    def __len__(self):
        return len(self.seq)

    def __getitem__(self, pos):
        return self.seq.__getitem__(pos)
        
