#!/usr/bin/env python3

import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import copy

from io import BytesIO
#import zipfile
import tarfile

sys.path.insert(0, '/home/')
import rpRanker


## Function that wraps the SBML files into a tar file to be passed to the next galaxy node
#
def writerpSBMLtar(rpsbml_paths, outTar):
    with tarfile.open(outTar, 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = BytesIO(data)
            info = tarfile.TarInfo(rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)


## Function that takes for input a tar of the collection of rpSBML files and reads them to memory
#
#
def readrpSBMLtar(inputTar):
    rpsbml_paths = {}
    tar = tarfile.open(inputTar)
    for member in tar.getmembers():
        rpsbml_paths[member.name] = rpRanker.rpSBML(member.name,libsbml.readSBMLFromString(tar.extractfile(member).read().decode("utf-8")))
    return rpsbml_paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Given an SBML model and the generated SBML heterologous pathway by RetroPath2.0, merge the two')
    parser.add_argument('-inModel', type=str)
    parser.add_argument('-inSBMLtar', type=str)
    parser.add_argument('-outSBMLtar', type=str)
    params = parser.parse_args()
    #sbml read the different mode
    rpsbml_paths = readrpSBMLtar(params.inSBMLtar)
    for rpsbml_name in rpsbml_paths:
        #read the input sbml model
        input_rpsbml = rpRanker.rpSBML('inputMergeModel')
        input_rpsbml.readSBML(params.inModel)
        rpsbml_paths[rpsbml_name].mergeModels(input_rpsbml.model)
        rpfba = rpRanker.rpFBA(input_rpsbml)
        rpfba.allObj()
        ##### pass FBA results to the original model ####
        groups = rpfba.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup('rp_pathway')
        for member in rp_pathway.getListOfMembers():
            reacFBA = rpfba.rpsbml.model.getReaction(member.getIdRef())
            reacIN = rpsbml_paths[rpsbml_name].model.getReaction(member.getIdRef())
            reacIN.setAnnotation(reacFBA.getAnnotation())
        #rpsbml_paths[rpsbml_name] = rpsbml_paths[rpsbml_name]
        input_rpsbml = None
    writerpSBMLtar(rpsbml_paths, params.outSBMLtar)
    exit(0)

