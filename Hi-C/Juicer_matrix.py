#! /usr/bin/env python3
import time
import math
import unittest
import re,sys,os
import numpy as np
import pandas as pd
from scipy import sparse



#step 1  juicer 提取矩阵
#start,end 设定
class trans_hic():

	def __init__(self):

		self.hic=""
		self.start=0
		self.end=100000
		self.juicer_tools=""
		self.outdir=""
		self.bin=10000
		self.genome='hg19'
		self.juicer_dump_mat=""
		self.hitc_matrix=""
		self.prefix="sam1"
	
	def extract_matrix(self): 
		chr=self.chr;start=self.start;end=self.end
		chrom=self.chr.replace('chr','')
		juicer_dump_mat="{}/{}_{}_{}_{}_dump.mat".format(self.outdir,self.prefix,chr,start,end)
		self.juicer_dump_mat=juicer_dump_mat
		region="{0}:{1}:{2} {0}:{1}:{2}".format(chrom,str(start),str(end))
		cmd="/opt/juicer/scripts/juicer_tools dump observed NONE {}  {}  BP {} {}".format(self.hic,region,self.bin,juicer_dump_mat)
		print(cmd)
		os.system(cmd)
	
	def reform_matrix(self): 
		#-----------HiTC matrix---------------------
		chr=self.chr;start=self.start;end=self.end;bin=self.bin;genome=self.genome
		self.hitc_matrix="{}/{}_{}_{}_{}_hitc.mat".format(self.outdir,self.prefix,chr,start,end)
		print('juicer dump matrix....')
		print(self.juicer_dump_mat)
		mat=pd.read_table(self.juicer_dump_mat,names=['frag1','frag2','contacts'])
		#print('matrix head.....')
		#print(mat.head())
		min=math.ceil(int(start)/bin)*bin
		max=int(int(end)/bin)*bin
		N=int(end/bin)-math.ceil(start/bin)+1
		#---------------------- add header --------------------------
		inddf=np.arange(N)
		headers_ref=[genome for x in inddf]
		bin_num_df=pd.Series(inddf).apply(lambda x : str(x))
		headers_ref=pd.Series(headers_ref)
		chromdf=pd.Series([chr for x in list(range(N))])
		startdf=pd.Series(inddf*bin+min)
		enddf=pd.Series((inddf+1)*bin+min)
		headers_suf=chromdf.str.cat(startdf.apply(lambda x :str(x)),sep=':')
		headers_suf=headers_suf.str.cat(enddf.apply(lambda x:str(x)),sep="-")
		headers=bin_num_df.str.cat([headers_ref,headers_suf],sep="|")
		headers=list(headers)

		mat['b1']=mat['frag1'].apply(lambda x: (x-min)/bin)
		mat['b2']=mat['frag2'].apply(lambda x: (x-min)/bin)
		counts=sparse.coo_matrix((mat['contacts'],(mat['b1'],mat['b2'])),shape=(N, N),dtype=float).toarray()
		diag_matrix=np.diag(np.diag(counts))
		counts=counts.T + counts
		counts=counts-diag_matrix-diag_matrix
		df=pd.DataFrame(counts)
		df.columns=headers
		df.index=headers
		#print('DataFrame.....')
		#print(df.head())
		df.to_csv(self.hitc_matrix,sep="\t")
		return df

	def z_score_norm(self):
		print('z-score normlizaion ....................')
		df=self.reform_matrix()
		print('befor zscore.......')
		print(df.head())
		dsc = pd.DataFrame(np.ravel(df)).describe(include=[np.number])
		df = (df - dsc.ix['mean',0])/dsc.ix['std',0]
		print('after zscore....')
		print(df.head())
		return df


class Test_trans(unittest.TestCase):
	def test_trans(self):
		trhic=trans_hic()
		trhic.outdir="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/Data_trans/fanxuning/B_TET_27/Annoroad_Lego_Browser/example/"
		trhic.hic="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/Data_trans/fanxuning/B_TET_27/batch_processiong/Data/hicFiles/D129/chr1/D129.chr1.hic"
		trhic.chr='chr1'
		trhic.start=62932570
		trhic.end=63564575
		trhic.extract_matrix()
		trhic.reform_matrix()


if __name__ == '__main__':
	unittest.main()

