#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:09:40 2020

@author: Joshua Sakai
"""

import sys
import cairocffi as cairo
import re
import Bio
import argparse
from Bio import SeqIO

def get_args():
    """argparse input files"""
    parser = argparse.ArgumentParser(description='file name to be used')
    parser.add_argument('-f', '--file', type=str, required=True, help='File path to fasta containing reads of interest')
    parser.add_argument('-m', '--motif', type=str, required=True, help='File path to fasta containing motifs separated by \n')
    parser.add_argument('-o', '--output', required=True, help='File name for the output SVG')
    return parser.parse_args()

#setting argparse values
args = get_args()
file = args.file
motif = args.motif
output= args.output

#setting global variables

#dictionary for converting the motifs as written to a regex
reg_ref = {"y": "[uct]", "g": "g", "c": "c", "t":"t", "a": "a", "u": "t", "n": "[uctag]", "Y": "[UCT]", "G": "G","T":"T", "C": "C", "A": "A", "U": "T", "N": "[UCTAG]"}

#Opening the motif file and sticking each motif as written into a list after converting it into a regex that will seach 1 nt at a time
mot=""
with open(motif, "r") as m:
    m_ref = []
    m_raw = []
    for line in m:
        x = line.split('\t')[0]
        y = x.strip()
        m_raw.append(y)
        for letter in y:
            mot = mot + reg_ref[letter].lower()
        m_ref.append("(?=("+mot+"))")
        mot = ""


max_len = 0
num_record=0
for record in SeqIO.parse('Figure_1.fasta', "fasta"):
    num_record+=1
    if len(record) > max_len:
        max_len = len(record)

#creating a surface to draw on
image = cairo.SVGSurface(str(output)+".svg", (max_len+100), (num_record*200))
#setting context as the output image
context = cairo.Context(image)

#setting a counter to move through a list of different RGB values
motif_color=0
#setting the starting x/y coordinates for the output image
point_x = 10
point_y = 75

#Setting color pallette
#Cyan, Blue, Violet, Purple-Red, Red, Light Gray, Dark Gray, Orange, Yellow, Green 
RGB = [[0,0.8,0.8],[0,0,0.8],[0.6,0,0.6],[0.8,0,0],[0.25,0.25,0.25],
       [0.5,0.5,0.5],[0.6,0,0.3],[0.8,0.4,0],[0.8,0,0.8],[0,0.8,0]]


for record in SeqIO.parse(file, "fasta"):   
    #for each record in the fasta draw a line where 1 pixel = 1 nt
    #for each reacord draw a box to represent the exon 
    header = str('>'+record.id)
    sequence = str(record.seq)
    back =str(record.seq)
    exon = re.findall(r'[A-Z]+', sequence)
    
    #setting aesthetics
    context.set_line_width(2)
    context.set_source_rgb(0, 0, 0)
    
    #header
    context.move_to(point_x, point_y)
    context.set_font_size(10)
    context.show_text(header)
    
    #drawing a line for the sequence
    context.move_to(point_x+50, point_y)
    context.line_to((len(sequence)+point_x+50), point_y)
    context.move_to(point_x, point_y-30)
    
    #drawing a rectangle for the exon
    exon_start = sequence.find(exon[0]) + 1 + point_x
    context.rectangle(exon_start, point_y-15, len(exon[0]), 30)
    context.stroke()

    for m in m_ref:
        #required to reset the sequence between motifs
        sequence = str(record.seq)
        first = True
        #tracking current position
        pos = 1
        #using the regex motifs, returns a list of all found instances of the motif in the sequence
        loc = re.findall(m, sequence.lower())

        #setting the color for all instances of each motif
        motif_color=m_ref.index(m)
        context.set_source_rgb(RGB[motif_color][0], RGB[motif_color][1], RGB[motif_color][2])
        context.stroke()
        
        #printing a legend value for a motif that appears 0 times in a sequence
        if len(loc)==0:
            context.set_font_size(10)
            context.move_to(point_x, point_y + (len(m_ref)*9) + (motif_color * 10))
            string= str("There are 0 occurances of "+ m_raw[motif_color]) + " in this record"
            context.show_text(string)
        
        for x in range(len(loc)):
            #now we know how many times the motif occurs we need to walk through the sequence and draw each occurance relative to its position in the sequence
            
            #find the first occurance of the motif
            motifs = sequence.lower().find(loc[x]) + pos

            #draw a line, width is 1pixel per 1nt
            context.set_line_width(1+ (len(m)/10))
            context.move_to(motifs+point_x+50 + (len(loc)/2), point_y+10)
            context.line_to(motifs+point_x+50 + (len(loc)/2), point_y-10)
            context.stroke()
            
            #cutting the sequence up until 1 nt past the motif while updating our position
            pos += sequence.lower().find(loc[x]) +1
            sequence = sequence[(sequence.lower().find(loc[x])+1):]    
            
            if first == True:
                #the first time a motif is found write out a legend showing sequence, color, and # of times the motif occured
                context.set_font_size(10)
                context.move_to(point_x, point_y + (len(m_ref)*9) + (motif_color * 10))
                string= str("There are " + str(len(loc)) + " occurances of "+ m_raw[motif_color]) + " in this record"
                context.show_text(string)
                #setting to False so this isn't repeated
                first = False

    #update starting y position for next record
    point_y += (len(m_ref)*25)


image.finish()





