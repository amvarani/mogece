#!/usr/bin/env python3
#
# MoGECE: extract coordinates for mobile genetic elements from outputs
# produced by MGE detection software
#
# Version 1.0 - March 19, 2016
#
# Copyright © 2016 Danillo Oliveira Alvarenga
#
# MoGECE is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# MoGECE is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more 
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with MoGECE. If not, see <http://www.gnu.org/licenses/agpl-3.0.html>.

import argparse

# Get command line arguments.
parser = argparse.ArgumentParser(
            description = "Mobile Genetic Elements Coordinates Extractor")

parser.add_argument("-v", "--version", action = "version",
                    version = "%(prog)s 1.0", help = "show version and exit")

parser.add_argument("-f", "--file", metavar = "File", required = True,
                    help = "prediction file")

group = parser.add_mutually_exclusive_group(required = True)
group.add_argument("-l", "--alienhunter", action = "store_true",
                   help = "using Alien_Hunter sco prediction file")
group.add_argument("-i", "--issaga", action = "store_true",
                   help = "using ISsaga csv prediction file")
group.add_argument("-m", "--minced", action = "store_true",
                   help = "using MinCED gff prediction file")
group.add_argument("-o", "--oasis", action = "store_true",
                   help = "using OASIS gff prediction file")
group.add_argument("-w", "--oligowords", action = "store_true",
                   help = "using OligoWords out prediction file")
group.add_argument("-p", "--phispy", action = "store_true",
                   help = "using PhiSpy tbl prediction file")
group.add_argument("-s", "--sniffer", action = "store_true",
                   help = "using SeqWordSniffer out prediction file")
group.add_argument("-r", "--virsorter", action = "store_true",
                   help = "using VirSorter fasta prediction file")

parser.add_argument("-a", "--artemis", action = "store_true",
                    help = "generate an Artemis gff file")
parser.add_argument("-g", "--gview", action = "store_true",
                    help = "generate a GView csv file")

args = parser.parse_args()

if not (args.gview or args.artemis):
    parser.error("Please indicate a visualization program")

# Create lists for MGE coordinates and scores
beginnings = []
ends = []
scores = []
isfamilies = []

# Get coordinates from an Alien_Hunter output file.
def alien_hunter():

    with open(args.file, "rt") as output:

        for line in output:

            line = line.split(' ')

            if "misc_feature" in line:
                beginning = line[-1].split("..")[0].strip()
                end = line[-1].split("..")[1].strip()
                add_coordinates(beginning, end)

            else:
                score = line[-1].split("=")[1].strip()
                score = str(round(float(score), 1))
                scores.append(score)

    return ("alienhunter", "genomic island", "Alien_Hunter")

# Get coordinates from an ISsaga output file.
def is_saga():

    with open(args.file, "rt") as output:

        IS = False
        ORF_L = 0
        ORF_R = 0
        family = 0

        for line in output:

            if "false positive" in line:
                continue

            line = line.split(',')

            if len(line) > ORF_R:

                if 'orf_L' and 'orf_R' in line:
                    ORF_L = line.index('orf_L')
                    ORF_R = line.index('orf_R')
                    family = line.index('family')

                elif ORF_L is not 0 and ORF_R is not 0:
                    beginning = line[ORF_L]
                    end = line[ORF_R]
                    add_coordinates(beginning, end)

                    scores.append('1')

                    isfamilies.append(line[family])

    return ("issaga", "transposase", "ISsaga")

# Get coordinates from a MinCED output file.
def min_ced():

    with open(args.file, "rt") as output:

        for line in output:

            if "##gff" in line:
                continue

            line = line.split("\t")

            beginning = line[line.index("CRISPR") + 1]
            end = line[line.index("CRISPR") + 2]
            add_coordinates(beginning, end)

            scores.append('1')

    return ("minced", "CRISPR", "MinCED")

# Get coordinates from an OASIS output file.
def oas_is():

    with open(args.file, "rt") as output:

        for line in output:

            line_items = list(line.split("\t"))

            beginning = line_items[3]
            end = line_items[4]
            add_coordinates(beginning, end)

            scores.append('1')

    return ("oasis", "transposase", "OASIS")

# Get coordinates from an OligoWords output file.
def oligo_words():

    with open(args.file, "rt") as output:

        island = False

        for line in output:

            if island:

                if not line.rstrip():
                    break

                line_items = list(line.split("\t"))

                beginning = line_items[0]
                end = line_items[1]
                add_coordinates(beginning, end)

                score = str(round(float(line_items[2]), 1))
                scores.append(score)

            if "n0_4mer:PS" in line:
                island = True

    return ("oligowords", "genomic island", "OligoWords")

# Get coordinates from a PhiSpy output file.
def phi_spy():

    with open(args.file, "rt") as output:

        for line in output:

            line_items = list(line.strip().split("\t"))

            beginning = line_items[1].split('_')[-2]
            end = line_items[1].split('_')[-1]
            add_coordinates(beginning, end)

            scores.append('1')

    return ("phispy", "prophage", "PhiSpy")

# Get coordinates from a SeqWord Gene Island Sniffer output file.
def seqword_sniffer():

    with open(args.file, "rt") as output:

        for line in output:

            if "<COORDINATES>" in line:

                line = line.split(' ')

                coordinates = line.index("<COORDINATES>") + 1               

                beginning = line[coordinates].split("-")[0]
                end = line[coordinates].split("-")[1]
                add_coordinates(beginning, end)

                score = str(round(float(line[-1].strip()), 1))
                scores.append(score)

    return ("sniffer", "genomic island", "SeqWord Sniffer")

# Get coordinates from a VirSorter output file.
def vir_sorter():

    with open(args.file, "rt") as output:

        for line in output:

            if '>' in line:

                line = line.split("_gene_")

                beginning = line[2].split("-")[1]
                end = line[2].split("-")[2]
                add_coordinates(beginning, end)

                scores.append('1')

    return ("virsorter", "prophage", "VirSorter")

# Add coordinates to lists making sure the ORF is correct.
def add_coordinates(beginning, end):

    if int(beginning) > int(end):
        beginning, end = end, beginning

    beginnings.append(beginning)
    ends.append(end)

# Generate a range file containing coordinates for visualization in GView.
def create_csv(feature):

    with open(feature[0] + ".csv", "wt") as csv:

        i = 0

        while i < len(scores):
            csv.write(beginnings[i] + ',' + ends[i] + ',' + scores[i] + '\n')
            i += 1

# Generate a feature table from coordinates for visualization in Artemis.
def create_ft(feature):

    with open(feature[0] + ".ft", "wt") as ft:

        i = 0

        while i < len(scores):

            if feature[1] == "transposase":
                ft.write("FT   repeat_region   " +
                         beginnings[i] + ".." + ends[i] + "\n")
            elif feature[1] == "prophage":
                ft.write("FT   misc_feature    " +
                         beginnings[i] + ".." + ends[i] + "\n")
            elif feature[1] == "CRISPR":
                ft.write("FT   misc_feature    " +
                         beginnings[i] + ".." + ends[i] + "\n")
            elif feature[1] == "genomic island":
                ft.write("FT   mobile_element  " +
                         beginnings[i] + ".." + ends[i] + "\n")

            if not isfamilies:
                ft.write("FT                   /note=" + feature[1] + "\n")
            else:
                if isfamilies[i]:
                    ft.write("FT                   /note=" + feature[1] +
                             " from the " + isfamilies[i].rstrip() +
                             " family\n")
                else:
                    ft.write("FT                   /note=" + feature[1] +
                             " from an unidentified family\n")

            ft.write("FT                   " +
                     "/annotation_source=" + feature[2] + "\n")

            if scores[i] is '1':
                i += 1
                continue
            else:
                ft.write("FT                   " +
                         "/score=" + scores[i] + "\n")
                i += 1

# Run functions according to corresponding arguments.
def main():

    if args.alienhunter:
        feature = alien_hunter()
    elif args.issaga:
        feature = is_saga()
    elif args.minced:
        feature = min_ced()
    elif args.oasis:
        feature = oas_is()
    elif args.phispy:
        feature = phi_spy()
    elif args.virsorter:
        feature = vir_sorter()
    elif args.sniffer:
        feature = seqword_sniffer()
    elif args.oligowords:
        feature = oligo_words()

    if args.artemis:
        create_ft(feature)
    if args.gview:
        create_csv(feature)

if __name__ == "__main__":
    main()
