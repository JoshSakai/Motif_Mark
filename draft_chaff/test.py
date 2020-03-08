# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 12:10:19 2020

@author: the5t
"""

import sys
import cairocffi as cairo
import re
import Bio
from Bio import SeqIO

m_ref = ['ygcy', 'GCAUG', 'catag', 'YYYYYYYYYY']

reg_ref = {"y": "[uct]", "g": "g", "c": "c", "t":"t", "a": "a", "u": "u", "n": "[uctag]", "Y": "[UCT]", "G": "G","T":"T", "C": "C", "A": "A", "U": "U", "N": "[UCTAG]"}

m_ref2 = []
mot = ""

for cha in m_ref:
    for letter in cha:
#        letter.replace(letter, reg_ref[letter])
        mot = mot + reg_ref[letter]
        print(mot)
    m_ref2.append(mot)
    mot = ""
#        m_ref= letter.replace(letter, reg_ref[letter])
        
#with open ('C:/bi625/Figure_1.fasta', 'r') as file:

#    for line in file:
#        if line.startswith(">"):
#            print(line)
#        else:
#            line.strip()
#            print(line)

image = cairo.PDFSurface("Figure_1.pdf", 750, 300)


def draw_motif(output, motif, move_x, move_y, line_x, line_y, c1, c2, c3):
    out = cairo.Context(output)
    out.set_line_width(len(motif))
    out.move_to(move_x + (len(motif)/2), move_y)
    out.line_to(line_x + (len(motif)/2), line_y)
    out.set_source_rgb(c1, c2, c3)
    out.stroke()


for record in SeqIO.parse('C:/bi625/Figure_1.fasta', "fasta"):
    context = cairo.Context(image)
    context.set_line_width(1)
    point_x = 10
    point_y = 50
    
    header = str('>'+record.id+'\n')
    sequence = str(record.seq)
    exon = re.findall(r'[A-Z]+', sequence)
    

    context.move_to(point_x, point_y)
    context.line_to((len(sequence)+point_x), point_y)
    context.set_font_size(7)
    context.move_to(point_x, point_y-13)
    context.show_text(header)
    
    exon_start = sequence.find(exon[0]) + 1 + point_x
    context.rectangle(exon_start, point_y-10, len(exon[0]), 20)
    
    for m in m_ref2:
        pos = 1
        loc = re.findall(m.lower(), sequence.lower())
        for x in range(len(loc)):
            motifs = sequence.lower().find(mot[x]) + loc
            print(motifs)
            
