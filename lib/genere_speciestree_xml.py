#
# Part of this  code is taken from https://groups.google.com/forum/#!topic/etetoolkit/cYTHHsL21KY
# written by Jaime Huerta-Cepas

"""Ecrit un  arbre au format XML .
Usage:
  genere_xml.py <treeFile>

Positional arguments:
  treeFile            hylogenic trees in newick format


"""
import sys
import re
import random
import os
from Bio import Phylo
from io import StringIO
# from cStringIO import StringIO
import time
from ete3 import Phyloxml, phyloxml

from lxml import etree
from xml.etree import ElementTree
from xml.dom import minidom

from lxml.etree import XMLParser, parse

from docopt import docopt
import sqlite3
import zlib
import base64



def createPhyloXML(newick):
    specoutfile= open("./instance/DGINN/results/species_list","w")
    handle = StringIO(newick)
    trees = Phylo.read(handle, 'newick')
    rd = str(random.randint(0,1000))
    Phylo.write([trees], 'tmpfile-'+rd+'.xml', 'phyloxml')
    file = open('tmpfile-'+rd+'.xml', 'r')
    text = file.read()
    file.close()
    os.remove('tmpfile-'+rd+'.xml')
    p = XMLParser(huge_tree=True)
    text = text.replace("phy:", "")
    text = re.sub("b'([^']*)'", "\\1", text)
    text = re.sub('branch_length_attr="[^"]+"', "", text)
    header = "<phyloxml>"

    text = re.sub('<phyloxml[^>]+>', header, text)
    text = text.replace('Phyloxml', 'phyloxml')
    tree = etree.fromstring(text,parser=p)
    # ajout du nom d'arbre
    treename = etree.Element("name")
    treename.text = "SPECIES_TREE"
    ins = tree.find('phylogeny')
    ins.append(treename)

    clade = tree.xpath("/phyloxml/phylogeny/clade")
    subtree = tree.xpath("/phyloxml")
    nbfeuille = 0
    famspecies = {}

    for element in clade[0].iter('clade'):
        enom=element.find('name')
        eleaves=element.find('clade')
        if ((enom is not None) and (eleaves is  None)) :
            nbfeuille = nbfeuille + 1
            cds = enom.text
            print(cds)
            specoutfile.write(cds+"\n")
            #sp = dico.get(cds)
            sp=cds;
            if (not  sp):
                print ("undefined species for "+ cds)
                sp = "undefined"
            famspecies[sp] = 1
            synLeft=""
            synRight = ""
            evrec = etree.Element("eventsRec")
            leaf = etree.Element("leaf")
            leaf.set('speciesLocation', sp)
            leaf.set('synteny', synLeft)
            #leaf.set('syntenyLeft', synLeft)
            #leaf.set('syntenyRight', synRight)
            evrec.append(leaf)
            element.append(evrec)
    print ("Number of leaves : ")
    print (nbfeuille)
    nbspecies = len(famspecies)
    print ("Number of species : ")
    print (nbspecies)

    treesize =  etree.Element("size")
    treesize.set('leaves',str(nbfeuille))
    treesize.set('species',str(nbspecies))
    e=subtree[0].find('phylogeny')
    e.append(treesize)
    text =  minidom.parseString(ElementTree.tostring(subtree[0])).toprettyxml()
    # remove blank lines
    cleantext = "\n".join([ll.rstrip() for ll in text.splitlines() if ll.strip()])
    return cleantext

