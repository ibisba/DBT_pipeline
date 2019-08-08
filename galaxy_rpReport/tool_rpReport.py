#!/usr/bin/env python3

import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import csv

from io import BytesIO
import tarfile

sys.path.insert(0, '/home/')
import rpRanker


## Function that takes for input a tar of the collection of rpSBML files and reads them to memory
#
#
def readrpSBMLtar(inputTar):
    rpsbml_paths = {}
    tar = tarfile.open(inputTar)
    for member in tar.getmembers():
        rpsbml_paths[member.name] = rpRanker.rpSBML(member.name,libsbml.readSBMLFromString(tar.extractfile(member).read().decode("utf-8")))
    return rpsbml_paths


def writeReport(rpsbml_name, rpsbml, csvfi, pathId='rp_pathway'):
    #loop through all the groups reactions
    groups = rpsbml.model.getPlugin('groups')
    rp_pathway = groups.getGroup(pathId)
    annot = rp_pathway.getAnnotation()
    path_annot = annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
    ### pathway
    #return all the names and select the ones that start with "fba_*"
    # --> duplicate the line if there are multiple
    toWrite = [rpsbml_name,
        'NaN',
        path_annot.getChild('dfG_prime_o').getAttrValue('value'),
        path_annot.getChild('dfG_prime_m').getAttrValue('value'),
        path_annot.getChild('dfG_uncert').getAttrValue('value'),
        path_annot.getChild('rule_id').getChild(0).toXMLString(),
        path_annot.getChild('smiles').getChild(0).toXMLString(),
        path_annot.getChild('rule_score').getAttrValue('value')]
    for i in range(path_annot.getNumChildren()):
        ent = path_annot.getChild(i)
        if ent.getName()[:4]=='fba_':
            toWrite.append(ent.getAttrValue('value'))
    csvfi.writerow(toWrite)
    ### reaction
    for member in rp_pathway.getListOfMembers():
        reaction = rpsbml.model.getReaction(member.getIdRef())
        annot = reaction.getAnnotation()
        react_annot = annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
        selen_annot = annot.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild('selenzyme')
        toWrite = [rpsbml_name,
            reaction.getId(),
            react_annot.getChild('dfG_prime_o').getAttrValue('value'),
            react_annot.getChild('dfG_prime_m').getAttrValue('value'),
            react_annot.getChild('dfG_uncert').getAttrValue('value'),
            react_annot.getChild('rule_id').getChild(0).toXMLString(),
            react_annot.getChild('smiles').getChild(0).toXMLString(),
            react_annot.getChild('rule_score').getAttrValue('value'),
            ';'.join([react_annot.getChild(i).getName() for i in range(react_annot.getNumChildren()) if react_annot.getChild(i).getName()[:4]=='fba_']),
            ';'.join([react_annot.getChild(i).getAttrValue('value') for i in range(react_annot.getNumChildren()) if react_annot.getChild(i).getName()[:4]=='fba_']),
            ';'.join([selen_annot.getChild(i).getName() for i in range(selen_annot.getNumChildren())]),
            ';'.join([selen_annot.getChild(i).getAttrValue('value') for i in range(selen_annot.getNumChildren())])]
        csvfi.writerow(toWrite)
        ### species
        #for pro in reaction.getListOfProducts():
        #for rea in reaction.getListOfReactants():
        #sel = bag_ibisba.getChild('selenzyme')


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Given an SBML, extract the reaction rules and pass them to Selenzyme REST service and write the results to the SBML')
    parser.add_argument('-inSBMLtar', type=str)
    parser.add_argument('-reportCSV', type=str)
    params = parser.parse_args()
    #sbml read the different mode
    rpsbml_paths = readrpSBMLtar(params.inSBMLtar)
    with open(params.reportCSV, 'w') as fi:
        csvfi = csv.writer(fi, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvfi.writerow(['Pathway', 'Reaction', 'dfG_prime_o', 'dfG_prime_m', 'dfG_uncert', 'Rule_id', 'Rule_SMILES', 'Rule_score', 'fba_obj', 'fba_flux', 'Selenzyme_uniprot', 'Selenzyme_score']) #need to add the FBA obj
        for rpsbml_name in rpsbml_paths:
            writeReport(rpsbml_name, rpsbml_paths[rpsbml_name], csvfi)
    exit(0)
