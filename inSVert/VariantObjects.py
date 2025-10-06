from abc import ABC, abstractmethod

'''
Structural Variant Objects used to build a VCF file
'''

class StructuralVariant(ABC):
    
    def __init__(self, chrom:str, pos:int, length:int, id:str):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.id = id 
        self.ref = "N"
        self.qual = "."
        self.filter = "PASS"

    @abstractmethod
    def get_alt(self) -> str:
        #return the alt field
        pass

    @abstractmethod
    def get_info(self) -> str:
        #return the info field
        # END = each SVTYPE processes END differently
        pass

    def format(self) -> str:
        #formats information into a VCF line
        alt = self.get_alt()
        info = self.get_info()
        return f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{alt}\t{self.qual}\t{self.filter}\t{info}\tGT\t1/1"
        

class Insertion(StructuralVariant):
#implement also sequence as an optional argument to pass tk get_alt

    def __init__(self, chrom, pos, length, id):
        super().__init__(chrom, pos, length, id)

    def get_alt(self):
        return "<INS>" #for now, I'll keep the symbolic alt
    
    def get_info(self):
        END = self.pos
        return f"SVTYPE=INS;SVLEN={self.length};END={END}"

ins1 = Insertion('ch1', 1000, 55, 'inSVert.INS.1')
print(ins1.format())

class Deletion(StructuralVariant):

    def __init__(self, chrom, pos, length, id):
        super().__init__(chrom, pos, length, id)

    def get_alt(self):
        return "<DEL>"

    def get_info(self):
        END = self.pos + self.length
        return f"SVTYPE=DEL;SVLEN={self.length};END={END}"
    
del1 = Deletion('chX',90, 800, 'inSVert.DEL.1')
print(del1.format())


class Inversion(StructuralVariant):

    def __init__(self, chrom, pos, length, id):
        super().__init__(chrom, pos, length, id)

    def get_alt(self):
        return "<INV>"

    def get_info(self):
        END = self.pos + self.length
        return f"SVTYPE=INV;SVLEN={self.length};END={END}"
    
inv1 = Inversion('chY',1000, 88, 'inSVert.INv.1')
print(inv1.format())



class Duplication(StructuralVariant):

    def __init__(self, chrom, pos, length, id, copy_number:int):
        super().__init__(chrom, pos, length, id)
        self.copy_number = copy_number

    def get_alt(self):
        return "<DUP>"
    
    def get_info(self):
        END = self.pos + self.length
        return f"SVTYPE=DUP;SVLEN={self.length};END={END}"
    
    def format(self) -> str:
        #overrides the format method of the parent class : adds details of DUPs 
        alt = self.get_alt()
        info = self.get_info()
        return f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{alt}\t{self.qual}\t{self.filter}\t{info}\tGT:CN\t1/1:{self.copy_number}"
    

dup1 = Duplication('chr1', 10000, 5000, 'inSVert.DUP.1', copy_number=4)
print(dup1.format())