import unittest
import vep_input_format_validator

class TestValidator(unittest.TestCase):

    #### Basic tests ####

    def test_valid_line(self):
        line = "5 140532 140532 T/C +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_valid_line - test a valid line")
        self.assertEqual(report, '')

    #### Strand ####

    def test_invalid_line_strand(self):
        line = "3 319781 319782 A/- f"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_strand - strand value/character not supported")
        self.assertEqual(report, "Strand is not in the expected format (+ or -)")

    #### Allele string ####

    def test_valid_line_allele_sv_char(self):
        line = "5 140532 140834 DEL +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_valid_line_allele_sv_char - test a valid SV entry")
        self.assertEqual(report, '')

    def test_invalid_line_allele_char(self):
        line = "5 140532 140532 T/0 +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_allele_char - for insertion, the start coordinate should be greater than the end coordinate")
        self.assertEqual(report, "Non supported characters in the allele 'T/0'")

    def test_invalid_line_allele_no_alt(self):
        line = "5 140532 140532 T +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_allele_no_alt - missing alternative allele")
        self.assertEqual(report, "Allele 'T' is not in the supported format 'REF/ALT' (e.g. A/C), except for structural variants")

    def test_invalid_line_allele_no_alt_2(self):
        line = "5 140532 140532 T/ +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_allele_no_alt_2 - missing alternative allele after '/'")
        self.assertEqual(report, "Allele 'T/' is not in the supported format 'REF/ALT' (e.g. A/C), except for structural variants")

    def test_invalid_line_allele_no_ref(self):
        line = "5 140532 140532 /C +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_allele_no_ref - missing reference allele")
        self.assertEqual(report, "Allele '/C' is not in the supported format 'REF/ALT' (e.g. A/C), except for structural variants")

    #### Coordinates ####

    def test_valid_line_insertion(self):
        line = "1 881907 881906 -/C +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_valid_line_insertion - Test a valid line with variant insertion")
        self.assertEqual(report, '')

    def test_invalid_line_insertion(self):
        line = "1 881906 881906 -/C +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_insertion - for insertion, the start coordinate should be greater than the end coordinate")
        self.assertEqual(report, "For insertion, the start coordinate should be greater than the end coordinate")

    def test_invalid_line_allele_length(self):
        line = "3 319781 319782 A/- +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_allele_length - allele length don't match the given coordinates")
        self.assertEqual(report, "Allele length don't match the given coordinates")

    #### Ordering ####

    def test_valid_line_ordered(self):
        vep_input_format_validator.previous_line = "5 130500 130500 G/C +"
        vep_input_format_validator.chrs_seen = { '5' : 1 }
        line = "5 140532 140532 T/C +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_valid_line_ordered - Variants ordered by location")
        self.assertEqual(report, '')

    def test_invalid_line_ordered(self):
        vep_input_format_validator.previous_line = "5 130500 130500 G/C +"
        vep_input_format_validator.chrs_seen = { '5' : 1 }
        line = "5 120500 120500 A/G +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_ordered - Variants not ordered by location")
        self.assertEqual(report, 'Not ordered: 5:120500 vs previous line 5:130500')

    def test_invalid_line_ordered_chr(self):
        vep_input_format_validator.previous_line = "5 130500 130500 G/C +"
        vep_input_format_validator.chrs_seen = { '5' : 1, '8' : 1 }
        line = "8 500255 500255 C/T +"
        report = vep_input_format_validator.check_line(line)
        print("\n> test_invalid_line_ordered_chr - Variants not ordered by chromosome block")
        self.assertEqual(report, "Not ordered: this entry on chromosome '8' is not ordered in the chromosome '8' block")

if __name__ == '__main__':
    unittest.main()

