# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 21:21:04 2022

@author: TingYi
"""

import argparse
import sinum_funcs
import pandas as pd

parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-f", type=str , default="example/dataset/ChuType_exampleGEM_log2.txt" , help="gene expression matrix (GEM) file")
parser.add_argument("-b", type=float , default=0.2 , help="box size")
parser.add_argument("-z", type=float , default=0.0 , help="z-score threshold")
parser.add_argument("-o", type=str , default="./example" , help="output directory")
parser.add_argument("-n", type=str , default="ChuType_example" , help="save name")

args = parser.parse_args()
GEM_file , bs , cf, outdir, savename = args.f, args.b, args.z, args.o, args.n

GEM=pd.read_csv(GEM_file,sep='\t',index_col=0)

sinum_funcs.sinumnet(GEM, outdir, savename, boxsize = bs, cutoff=cf)

