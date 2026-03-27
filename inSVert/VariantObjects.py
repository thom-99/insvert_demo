from abc import ABC, abstractmethod

'''
Structural Variant Objects used to build a VCF file
'''

class StructuralVariant(ABC):
    
    def __init__(self, chrom:str, pos:int, length:int, id:str, genotype:str):
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.id = id 
        self.genotype = genotype
        self.ref = "N"
        self.qual = "."
        self.filter = "PASS"

    @abstractmethod
    def get_alt(self) -> str:
        #return the alt field
        pass

    @abstractmethod
    def get_end(self) -> str:
        # each svtype processes END differently
        pass

    @abstractmethod
    def get_info(self) -> str:
        #return the info field
        pass

    def format(self) -> str:
        #formats information into a VCF line
        alt = self.get_alt()
        info = self.get_info()
        return f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{alt}\t{self.qual}\t{self.filter}\t{info}\tGT\t{self.genotype}"
        

class Insertion(StructuralVariant):
#implement also sequence as an optional argument to pass tk get_alt

    def __init__(self, chrom, pos, length, id, genotype):
        super().__init__(chrom, pos, length, id, genotype)

    def get_alt(self):
        return "<INS>" #for now, I'll keep the symbolic alt
    
    def get_end(self):
        return self.pos
    
    def get_info(self):
        END = self.get_end()
        return f"SVTYPE=INS;SVLEN={self.length};END={END}"



class Deletion(StructuralVariant):

    def __init__(self, chrom, pos, length, id, genotype):
        if length>= 0:
            raise ValueError(f"Deletion length must be negative, got {length} instead")
        super().__init__(chrom, pos, length, id, genotype)

    def get_alt(self):
        return "<DEL>"
    
    def get_end(self):
        return self.pos + abs(self.length)

    def get_info(self):
        END = self.get_end()
        return f"SVTYPE=DEL;SVLEN={self.length};END={END}"
    


class Inversion(StructuralVariant):

    def __init__(self, chrom, pos, length, id, genotype):
        super().__init__(chrom, pos, length, id, genotype)

    def get_alt(self):
        return "<INV>"

    def get_end(self):
        return self.pos + self.length

    def get_info(self):
        END = self.get_end()
        return f"SVTYPE=INV;SVLEN={self.length};END={END}"
    



class Duplication(StructuralVariant):

    def __init__(self, chrom, pos, length, id, genotype, copy_number:int):
        super().__init__(chrom, pos, length, id, genotype)
        self.copy_number = copy_number

    def get_alt(self):
        return "<DUP>"
    
    def get_end(self):
        return self.pos + self.length

    def get_info(self):
        END = self.get_end()
        return f"SVTYPE=DUP;SVLEN={self.length};END={END}"
    
    def format(self) -> str:
        #overrides the format method of the parent class : adds details of DUPs 
        alt = self.get_alt()
        info = self.get_info()
        return f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{alt}\t{self.qual}\t{self.filter}\t{info}\tGT:CN\t{self.genotype}:{self.copy_number}"
    


class Breakend(StructuralVariant):
    """
    Represents a single BND record in VCF 4.2.
    A full translocation is composed of 4 Breakend objects.
    """
    def __init__(self, chrom, pos, id, genotype, mate_id, event_id, alt_string, role):
        # Length is typically 0 or 1 for BNDs as they represent a single point junction
        super().__init__(chrom, pos, 0, id, genotype)
        self.mate_id = mate_id
        self.event_id = event_id
        self.alt_string = alt_string # The bracketed ALT string (e.g., N[chr2:5000[)
        self.role = role             # SOURCE or DESTINATION

    def get_alt(self):
        """Returns the specific bracketed ALT string for this breakend."""
        return self.alt_string

    def get_end(self):
        """For BNDs, the END is the same as the POS."""
        return self.pos

    def get_info(self):
        """
        Returns the INFO field including BND-specific tags.
        """
        # We include MATEID to link to the partner BND and EVENT to group all 4.
        # TRA_ROLE helps the simulation engine know if it's cutting or pasting.
        return (f"SVTYPE=BND;MATEID={self.mate_id};"
                f"EVENT={self.event_id};TRA_ROLE={self.role}")

    def format(self):
        """Overrides format to ensure SVLEN is not typically included for BNDs."""
        alt = self.get_alt()
        info = self.get_info()
        return (f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{alt}\t"
                f"{self.qual}\t{self.filter}\t{info}\tGT\t{self.genotype}")



'''
bnd1 = Breakend('chr1', 10, 'inSVert.BND.1.1', '0/1', 'inSVert.BND.1.3', 'TRA1', 'N[chr1:10[' ,'SOURCE')
print(bnd1.format())
'''