"""
   Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
   Copyright [2016-2021] EMBL-European Bioinformatics Institute
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import sys, os, re

## Methods ##

def check_line(line):
    """ Check the VEP input line format. Return error message if the format is not as expected """
    global col_sep, previous_line

    line_parts = re.split(col_sep, line)

    # Check column numbers
    if (len(line_parts) < 4):
        return "Missing column(s)"
    elif (len(line_parts) > 6):
        return "Too many column(s) ("+str(len(line_parts))+" columns in total)"

    chr    = line_parts[0]
    start  = int(line_parts[1])
    end    = int(line_parts[2])
    allele = line_parts[3]
    strand = line_parts[4]

    # Check line order
    previous_line_report = compare_with_previous_line(chr,start,end)
    if (previous_line_report):
        return previous_line_report

    # Check strand
    strand_report = check_strand(strand)
    if (strand_report):
        return strand_report

    # Check alleles
    allele_report = check_alleles(allele,start,end)
    if (allele_report):
        return allele_report

    return ''


def check_strand(strand):
    """ Check the strand format. Return error message if the format is not as expected. """

    if (strand != '-' and strand != '+'):
        return "Strand is not in the expected format (+ or -)"


def check_alleles(allele,start,end):
    """ Check the allele format, in several ways. Return error message if the format is not as expected. """

    # Check allele characters
    if (not(re.search('^[AT\/GC-]+$',allele)) and not(re.search('^(INS|DEL|TDUP|DUP)$',allele))):
        return "Non supported characters in the allele '"+allele+"'"
    # Check the format REF/ALT, except for structural variants
    elif (not(re.search('^(INS|DEL|TDUP|DUP)$',allele)) and not(re.search('[ATCG-]\/[ATCG-]',allele))):
        return "Allele '"+allele+"' is not in the supported format 'REF/ALT' (e.g. A/C), except for structural variants"

    alleles = allele.split('/')
    coord_length = end - start + 1

    # Insertion
    if ((alleles[0] == '-' or allele == 'INS') and coord_length != 0):
        return "For insertion, the start coordinate should be greater than the end coordinate"
    elif (len(alleles[0]) != coord_length and coord_length != 0 and not(re.search('^(INS|DEL|TDUP|DUP)$',alleles[0]))):
        return "Allele length don't match the given coordinates"


def compare_with_previous_line(chr,start,end):
    """ Compare the current line with the previous line to check the variants ordering (chromosome and position).
        Return error message if the variants are not ordered. """
    global col_sep, previous_line, chrs_seen

    if (previous_line != ''):

        min_coord = start > end and end or start

        previous_line_parts = re.split(col_sep, previous_line)
        previous_chr = previous_line_parts[0]
        previous_start = previous_line_parts[1] > previous_line_parts[2] and previous_line_parts[2] or previous_line_parts[1]
        if (chr == previous_chr and int(previous_start) > min_coord):
            return "Not ordered: "+chr+":"+str(min_coord)+" vs previous line "+previous_chr+":"+previous_start
        elif (chr != previous_chr and chrs_seen.get(chr)):
            return "Not ordered: this entry on chromosome '"+chr+"' is not ordered in the chromosome '"+chr+"' block"
        else:
            chrs_seen[chr] = 1


def main():
    """ Main method fetching the VEP input file, reading it line by line and checking the format of each of them. """
    global msg, previous_line

    # File path missing in command line
    if len(sys.argv) < 2:
      sys.exit("ERROR: missing options!\n\nPlease, use the following options:\n"+msg)

    vep_input_filepath = sys.argv[1]

    # File doesn't exist
    if not(os.path.isfile(vep_input_filepath)):
        sys.exit("ERROR: file '"+vep_input_filepath+"' doesn't exist")

    error_count = 0

    # Open input file
    with open (vep_input_filepath, "r") as fileHandler:
        # Read/parse each line of the input file
        line_number = 1
        for line in fileHandler:
            line_content = line.strip();
            check_status = check_line(line_content)
            # Formatting issue found in the current line
            if (check_status != ''):
                print("Line "+str(line_number)+": "+check_status+" | Line content: \""+line_content+"\"")
                error_count = error_count + 1
            previous_line = line_content
            line_number = line_number + 1

    # Print final summary report
    if (error_count == 0):
        print("\n=> Your VEP input file is valid!")
    else:
        print("\n=> You have "+str(error_count)+" lines not valid")


## Global variables ##

msg = '''  python3 vep_input_format_validator.py <path_to_vep_input_file>

This script will check:
  - The number of columns on each lines
  - The alleles: make sure it's limited to A,T,G,C,-,INS,DEL,TDUP,DUP
  - The allele string format is 'REF/ALT' (except for structural variants)
  - The strand ('+' or '-')
  - The coordinates matches the length of the reference alleles (except for structural variants)
  - In case of insertion, it will check that the start is greater than the end
  - The ordering of the variants by locations
  - All the variants mapped to a chromosome are in the same block

VEP input format definition can be found here: https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default
'''

chrs_seen = {}

col_sep = '\s+'

previous_line = ''

## Main ##
if __name__ == "__main__":
    main()

