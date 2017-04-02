#! /usr/bin/env python2

import vcf
import vcf.utils
import hgvs

__author__ = 'Priska Lang'


class Assignment3:
    """
        Provides code to calculate and print variants information of the files AmpliseqExome.20141120.NA24385.vcf,
            AmpliseqExome.20141120.NA24143.vcf and AmpliseqExome.20141120.NA24149.vcf.

            This class parses the three vcf files and calculates the following properties unsing
            the PyVCF module:
            - total number of variants of mother
            - total number of variants of father
            - number of variants shared by father and son
            - number of variants shared by mother and son
            - number of variants shared by mother, father and son
            - variants shared by mother, father  and son
            - convert the first 100 variants of son into HGVS
            ! The last method doesn't work!
            For more information see: https://pyvcf.readthedocs.io/en/latest/ and:
            https://hgvs.readthedocs.io/en/master/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript
            It provides methods for getting the listed properties above and "print_summary" for printing them.
            """
    def __init__(self):
        """
            The constructor method parses the three vcf files and validates if PyVCF and hgvs are installed and
            """
        self.vcfSon = vcf.Reader(open("AmpliseqExome.20141120.NA24385.vcf", "r"))
        self.vcfMother = vcf.Reader(open("AmpliseqExome.20141120.NA24143.vcf", "r"))
        self.vcfFather = vcf.Reader(open("AmpliseqExome.20141120.NA24149.vcf", "r"))
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        ## Check if hgvs is installed
        print("HGVS version: %s" % hgvs.__version__)

    def get_total_number_of_variants_mother(self):
        """
            Prints the total number of identified variants in the mother
            """
        numVariantsMother = 0
        for record in self.vcfMother:
            numVariantsMother += 1
        print("total number variants mother:")
        print(numVariantsMother)

    def get_total_number_of_variants_father(self):
        """
            Prints the total number of identified variants in the father
            """
        numVariantsFather = 0
        for record in self.vcfFather:
            numVariantsFather += 1
        print("\ntotal number of variants father:")
        print(numVariantsFather)

    def get_variants_shared_by_father_and_son(self):
        """
            Prints the number of identified variants shared by father and son. Uses vcf.utils method walk_together.
            For more information see: http://pyvcf.readthedocs.io/en/latest/_modules/vcf/utils.html#walk_together
            """
        vcfSon = vcf.Reader(open("AmpliseqExome.20141120.NA24385.vcf", "r"))
        vcfFather = vcf.Reader(open("AmpliseqExome.20141120.NA24149.vcf", "r"))
        sharedVariantsFatherSon = 0
        shared_variants = vcf.utils.walk_together(vcfFather, vcfSon)
        for record in shared_variants:
            if not record[0] is None and not record[1] is None:
                sharedVariantsFatherSon += 1
        print("\nnumber of variants shared by father and son:")
        print(sharedVariantsFatherSon)

    def get_variants_shared_by_mother_and_son(self):
        """
            Prints the number of identified variants shared by mother and son. Uses vcf.utils method walk_together.
            For more information see: http://pyvcf.readthedocs.io/en/latest/_modules/vcf/utils.html#walk_together
            """
        vcfSon = vcf.Reader(open("AmpliseqExome.20141120.NA24385.vcf", "r"))
        vcfMother = vcf.Reader(open("AmpliseqExome.20141120.NA24143.vcf", "r"))
        sharedVariantsMotherSon = 0
        shared_variants = vcf.utils.walk_together(vcfMother, vcfSon)
        for record in shared_variants:
            if not record[0] is None and not record[1] is None:
                sharedVariantsMotherSon += 1
        print("\nnumber of variants shared by mother and son:")
        print(sharedVariantsMotherSon)

    def get_variants_shared_by_trio(self):
        """
            Prints the number of identified variants shared by father, mother and son. Uses vcf.utils method walk_together.
            For more information see: http://pyvcf.readthedocs.io/en/latest/_modules/vcf/utils.html#walk_together
            """
        vcfSon = vcf.Reader(open("AmpliseqExome.20141120.NA24385.vcf", "r"))
        vcfMother = vcf.Reader(open("AmpliseqExome.20141120.NA24143.vcf", "r"))
        vcfFather = vcf.Reader(open("AmpliseqExome.20141120.NA24149.vcf", "r"))
        sharedVariantsTrio = 0
        shared_variants = vcf.utils.walk_together(vcfMother, vcfFather, vcfSon)
        for record in shared_variants:
            if not record[0] is None and not record[1] is None and not record[2] is None:
                sharedVariantsTrio += 1
        print("\nnumber of variants shared by mother, father and son:")
        print(sharedVariantsTrio)

    def merge_mother_father_son_into_one_vcf(self):
        """
            Creates one vcf file, named mergedAmlpseqExome.vcf containing all variants of the trio. Uses vcf.utils
            method walk_together.
            For more information see: http://pyvcf.readthedocs.io/en/latest/_modules/vcf/utils.html#walk_together
            """
        vcfSon = vcf.Reader(open("AmpliseqExome.20141120.NA24385.vcf", "r"))
        vcfMother = vcf.Reader(open("AmpliseqExome.20141120.NA24143.vcf", "r"))
        vcfFather = vcf.Reader(open("AmpliseqExome.20141120.NA24149.vcf", "r"))
        vcf_writer = vcf.Writer(open("mergedAmpliseqExome.vcf","w"), vcfSon)
        shared_variants = vcf.utils.walk_together(vcfMother, vcfFather, vcfSon)
        for record in shared_variants:
            if not record[0] is None and not record[1] is None and not record[2] is None:
                vcf_writer.write_record(record[0])
        print("\nmerged vcf file saved as mergedAmpliseqExome.vcf")

    def convert_first_variants_of_son_into_HGVS(self):
        '''
        ! this method is not implemented!
        Converts the first 100 variants identified in the son into the corresponding transcript HGVS.
        Each variant should be mapped to all corresponding transcripts. Pointer:
        https://hgvs.readthedocs.io/en/master/examples/manuscript-example.html#project-genomic-variant-to-a-new-transcript
        Five Types of variants in HGVS: http://varnomen.hgvs.org/bg-material/simple/
        HGVS format: reference:description
            eg: NM_004006.2:c.4375C>T
            In the example NM_004006.2 is the reference sequence and c.4375C>T is the description of a variant.
        Writes the output to hgvs_son.txt
        '''
        outfile = open("hgvs_son.txt","w")
        vcfSon = vcf.Reader(open("AmpliseqExome.20141120.NA24385.vcf", "r"))
        for i in range(100):
            record = vcfSon.next()
            if record.is_snp:                   # lt. HGVS: substitution -> >
                None
                # TODO: implement
            elif record.var_subtype is "del":   # lt. HGVS: deletion -> del
                None
                # TODO: implement
            elif record.var_subtype is "ins":   # lt. HGVS: insertion -> ins
                None
                # TODO: implement
            elif record.var_subtype is "dup":   # lt. HGVS: duplication -> dup
                None
                # TODO implement
            else:                               # lt. HGVS: deletion/insertion (indel) -> delins
                None
                # TODO implement

    def print_summary(self):
        """
            print_summary calls all methods and prints the results calculated within the methods
            self.get_total_number_of_variants_mother, self.get_total_number_of_variants_father,
            self.get_variants_shared_by_father_and_son, self.get_variants_shared_by_mother_and_son and
            self.get_variants_shared_by_trio.

                    :Example:

                    For the given input file the following is printed on the console:

                :Example:
                PyVCF version: 0.6.8
                HGVS version: 0.4.13
                total number variants mother:
                38693
                total number of variants father:
                38641
                number of variants shared by father and son:
                30142
                number of variants shared by mother and son:
                30216
                number of variants shared by mother, father and son:
                22533
                """
        self.get_total_number_of_variants_mother()
        self.get_total_number_of_variants_father()
        self.get_variants_shared_by_father_and_son()
        self.get_variants_shared_by_mother_and_son()
        self.get_variants_shared_by_trio()
        self.merge_mother_father_son_into_one_vcf()
        self.convert_first_variants_of_son_into_HGVS()


if __name__ == '__main__':
    print("Assignment 3")
    assignment1 = Assignment3()
assignment1.print_summary()