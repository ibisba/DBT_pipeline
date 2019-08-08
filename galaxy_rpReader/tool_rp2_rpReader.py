#!/usr/bin/env python3

import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import csv
import json

from io import BytesIO
#import zipfile
import tarfile

sys.path.insert(0, '/home/')
import rpRanker

## Function that wraps the SBML files into a zip file to be passed to the next galaxy node
#
def writerpSBMLtar(rpsbml_paths, outTar):
    with tarfile.open(outTar, 'w:xz') as tf:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            fiOut = BytesIO(data)
            info = tarfile.TarInfo(rpsbml_name)
            info.size = len(data)
            tf.addfile(tarinfo=info, fileobj=fiOut)


'''
def writerpSBMLzip(rpsbml_paths, outZip):
    zip_buffer = io.BytesIO()
    #with zipfile.ZipFile(zip_buffer, mode="a", compression=zipfile.ZIP_BZIP2) as zip_file:
    #with zipfile.ZipFile(zip_buffer, mode="a") as zip_file:
    with zipfile.ZipFile(zip_buffer, "a", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_file:
        for rpsbml_name in rpsbml_paths:
            data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
            data = io.BytesIO(bytes(data))
            zip_file.writestr(rpsbml_name, data.getvalue())
    with open(outZip, 'wb') as f:
        f.write(zip_buffer.getvalue())
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to read the output of rp2paths (or JSON) into SBML files')
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2paths_outPaths', type=str)
    parser.add_argument('-rp2paths_scope', type=str)
    parser.add_argument('-outSBMLtar', type=str)
    parser.add_argument('-maxRuleIds', type=int)
    params = parser.parse_args()

    rpreader = rpRanker.rpReader()
    
    rpreader.compounds(params.rp2paths_compounds)
    rpreader.transformation(params.rp2paths_scope)
    rpreader.outPaths(params.rp2paths_outPaths, params.maxRuleIds)
    rpreader.pathsToSBML()
    
    #rpcofactors = rpRanker.rpCofactors(rpreader)
    #rpcofactors.addCofactors()
    #rpcofactors.pathsToSBML()

    #package the rpsbml's for the next pathway
    writerpSBMLtar(rpreader.sbml_paths, params.outSBMLtar)
    #writerpSBMLtar(rpcofactors.sbml_paths, params.outSBMLtar)
    #writerpSBMLzip(rpcofactors.sbml_paths, params.outSBMLzip)
    exit(0)
