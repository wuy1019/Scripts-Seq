# conding: utf-8
__author__ = "wuy"

import re, os
import string
from collections import defaultdict
from optparse import OptionParser

def get_para():
    usage = """python %prog -g <gff file> -f <genome> -l <length> [options]
    --------------------------------------------------------------------
    Veision: 0.3, support NCBI gff file
    Author:  wuy1019@live.com
    Notice:  fasta file title must be equal column one of gff file
    --------------------------------------------------------------------"""
    parser = OptionParser(usage = usage)
    parser.add_option("-g", "--gff", action = "store",
                    # dest="gff",
                    help = "annotation file [gff format]")
    parser.add_option("-f", "--fasta", action = "store",
                    # dest="fasta",
                      help = "genome file [fasta format]")
    
    parser.add_option("-l","--length",action = "store",
                      # dest = "len",
                      help = "gene up length [bp]")
    parser.add_option("-c", "--coerce", action = "store",
                      # dest="coerce",
                      help = "ignore up gene, coerce obtain gene up region(yes or no)|[default: no]")

    (options, args) = parser.parse_args()
    return options, args
options, args = get_para()

Gff_f   = options.gff
Fasta_f = options.fasta
Length  = options.length
Coerce  = options.coerce

def Ddict():
    return defaultdict(dict)

def GffInfo(Gff_file):
    """

    :rtype : object
    :param:  Gff_file: gff format file
    :return: gff[contig][n] = [title, start, end, direction]
    """
    global i
    gff = Ddict()
    try:
        file_gff = open(Gff_file, "r")
        try:
            contig_set = {}
            for line in file_gff:
                line = line.strip()
                if re.match("^#", line):
                    continue
                line_l = re.split("\t", line)
                if line_l[2] == "CDS" or "RNA" in line_l[2]:
                    if line_l[0] not in contig_set.keys():
                        i = 1
                        contig_set[line_l[0]] = ""
                    else:
                        i += 1
                    cache_dict = dict([(item.split("=")[0], item.split("=")[1])
                                       for item in line_l[-1].split(";")])
                    # rele = re.compile(r"^ID=(.+);Name=")
                    gff[line_l[0]][i] = [cache_dict["ID"], line_l[3],
                                        line_l[4], line_l[6]]
                                         
        except:
            pass
        finally:
            file_gff.close()
    except IOError:
        print "cant't open gff file"
        pass
    return gff

def FastaInfo(fasta_file):
    "store fasta file"
    fasta = {}
    title = ""
    seq = ""
    try:
        file_fasta = open(fasta_file, "r")
        if re.match("^>", file_fasta.readline()):
            file_fasta.seek(0)
            for line in file_fasta:
                line = line.strip()
                if  re.match("^>", line) and seq:
                    fasta[title] = seq
                    title = line[1:]
                    seq   = ""
                elif re.match("^>", line):
                    title = line[1:]
                else:
                    seq += line
            fasta[title] = seq
        else:
            print "file is not fasta fortmat"
        file_fasta.close()
    except IOError:
        print "cant't open fasta file"
        pass
    return fasta

def GffEnd(Gff_file, fasta_file):
    """add up gene end site;
    gff_end[contig][gene_id] = ["+", up_gene_end-1,   start-1, end]
    gff_end[contig][gene_id] = ["-", down_gene_start-1, end, start-1]
    """
    gff_raw = GffInfo(Gff_file)
    gff_end = Ddict()
    for key1 in gff_raw.keys():
        gene_num = max(gff_raw[key1].keys())
        for key2 in sorted(gff_raw[key1].keys()):
            if  gff_raw[key1][key2][3] == "+":
                if  key2 == 1:
                    gff_end[key1][gff_raw[key1][key2][0]] = ["+", 0,
                    int(gff_raw[key1][key2][1]) - 1, int(gff_raw[key1][key2][2])]
                else:
                    gff_end[key1][gff_raw[key1][key2][0]] = ["+", 
                    int(gff_raw[key1][key2 - 1][2]) - 1, int(gff_raw[key1][key2][1]) - 1,
                    int(gff_raw[key1][key2][2])]
                    
            else:
                if  key2 == gene_num:
                    gff_end[key1][gff_raw[key1][key2][0]] = ["-", 
                    len(FastaInfo(fasta_file)[key1]), int(gff_raw[key1][key2][2]),
                    int(gff_raw[key1][key2][1]) - 1]
                else:
                    gff_end[key1][gff_raw[key1][key2][0]] = ["-",
                    int(gff_raw[key1][key2 + 1][1]) - 1, int(gff_raw[key1][key2][2]),
                    int(gff_raw[key1][key2][1]) - 1]
    return gff_end

def Complement(char):
    "reverse compliment gene seq"
    char = re.sub("\s+", "", char)
    libary = string.maketrans('atcgATCG','tagcTAGC')
    end = char[::-1].translate(libary)
    return end

def GeneUp(seq, direction, up_gene, start, end, length = 500, coerce = "no"):
    "obtain gene up region"
    # gff = GffEnd(Gff_file, fasta_file)
    if coerce == "yes":
        if direction == "+":
            if  start <= length:
                up_char, gene_seq = seq[ : start], seq[start : end]
            else:
                up_char, gene_seq = seq[ start - length: start], seq[start : end]
        else:
            if  len(seq) - start <= length:
                up_char, gene_seq = \
                Complement(seq[start : len(seq)]), Complement(seq[end : start])
            else:
                up_char, gene_seq = \
                Complement(seq[start : start + length]), Complement(seq[end : start])

    else:
        if direction == "+":
            if  start - up_gene <= length:
                up_char, gene_seq = seq[up_gene : start], seq[start : end]
            else:
                up_char, gene_seq = seq[start - length : start], seq[start : end]
        else:
            if  up_gene - start <= length:
                up_char, gene_seq = Complement(seq[start : up_gene]), \
                Complement(seq[end : start])
            else:
                up_char, gene_seq = Complement(seq[start : start + length]),\
                Complement(seq[end : start])
    return up_char, gene_seq

def EndInfo(Gff_file, fasta_file, length = 500, coerce = "no"):
    gff_h = GffEnd(Gff_file, fasta_file)
    gene_data = {}
    for key1 in gff_h.keys():
        seq = FastaInfo(fasta_file)[key1]
        for key2 in gff_h[key1].keys():
            gene_id   = key2
            direction = gff_h[key1][key2][0]
            up_gene   = gff_h[key1][key2][1]
            start     = gff_h[key1][key2][2]
            end       = gff_h[key1][key2][3]
            up_char   = GeneUp(seq, direction, up_gene, start, end, length, coerce)[0]
            gene_seq  = GeneUp(seq, direction, up_gene, start, end, length, coerce)[1]
            gene_data[gene_id] = [up_char, gene_seq]
    return  gene_data

def sorted_hash(hash):
    "sorted for hash key, prefix alphabet at first then digit tag"
    key_list = hash.keys()
    key_list.sort(key=lambda x : (re.split("(?<=\D)\d+$", x)[0],
                              int(re.split("\D(?=\d+$)", x)[-1])))
    return key_list

try:
    END_gene = open( os.path.splitext(Fasta_f)[0]+"_genes.fasta","w" )
    END_up = open( os.path.splitext(Fasta_f)[0]+"_up"+Length+".fasta","w")
    Data = EndInfo(Gff_f, Fasta_f, length=int(Length), coerce=Coerce)
    for id in sorted_hash(Data):
        if len(Data[id][0]) % 60 != 0:
            END_up.write(  ">%s\n%s\n"%(id, re.sub("(\w{60})", "\\1\n", Data[id][0])))
        else:
            END_up.write(  ">%s\n%s"%(id, re.sub("(\w{60})", "\\1\n", Data[id][0])))
        if len(Data[id][1]) % 60 != 0:
            END_gene.write(">%s\n%s\n"%(id, re.sub("(\w{60})", "\\1\n", Data[id][1])))
        else:
            END_gene.write(">%s\n%s"%(id, re.sub("(\w{60})", "\\1\n", Data[id][1])))
    END_up.close()
    END_gene.close()
except:
    os.system("python %s -h"%os.path.basename(__file__))
    pass

