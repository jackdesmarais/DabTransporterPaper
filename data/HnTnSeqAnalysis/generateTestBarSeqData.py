#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
from optparse import OptionParser
# from scipy import stats
from scipy.stats import skewnorm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from operator import add
import time
from os import path
from os import makedirs
# skewnorm=stats.norm
random_state=int(time.gmtime().tm_sec)

def makeFit_logratiosValues(df,numGenes,fractionImportant,noise):
	normNum=int(np.ceil(numGenes*(1-fractionImportant)))
	phenNum = numGenes-normNum
	c1Vals = skewnorm.rvs(-2, scale=0.99, size=normNum,random_state=random_state)
	c1PhenoVals = skewnorm.rvs(0, scale=noise, size=phenNum,random_state=random_state)
	c2Vals = list(map(add, c1Vals, skewnorm.rvs(0, scale=noise*2, size=normNum,random_state=random_state)))
	c2PhenoVals = skewnorm.rvs(-20, scale=4, size=phenNum,random_state=random_state)
	df['setAS2 c1']=np.append(c1Vals,c1PhenoVals)
	df['setAS3 c2']=np.append(c2Vals,c2PhenoVals)
	return(df)

def makeGeneName(num,digs):
	return 'gene_'+str(num).zfill(digs)

def makefit_logratiosdf(numGenes):
	digs = len(str(numGenes))
	genelist = list(map(makeGeneName,list(range(0,numGenes+1)),[digs]*numGenes))
	df = pd.DataFrame(columns=['locusId','sysName','desc'])
	df.locusId,df.sysName,df.desc = genelist,genelist,genelist
	return(df)

def makeFit_logratios(options):
	df = makefit_logratiosdf(options.numGenes)
	df = makeFit_logratiosValues(df,options.numGenes,options.fractionImportant,options.noise)
	df.to_csv(path.join(options.outDir, 'fit_logratios.tab'),sep='\t',index=False)
	fig = plt.figure(figsize=(10,10))
	plt.scatter(df['setAS2 c1'], df['setAS3 c2'], c='g', s=100)
	plt.savefig(path.join(options.outDir, 'fitnessGraph'))

def main(options, args):
	makeFit_logratios(options)


if __name__ == '__main__':
   parser = OptionParser()
   parser.add_option("-o", "--outDir", dest="outDir", default='./testTnSeqData', help="directory for output")
   parser.add_option("-g", "--numGenes", dest="numGenes",type='int', default=2000, help="Number of genes in the sample")
   parser.add_option("-f", "--fractionImportant", dest="fractionImportant",type='float', default=0.01, help="fraction of genes with a phenotype")
   parser.add_option("-n", "--noise", dest="noise",type='float', default=0.12, help="how much noise to add to the data")

   (options, args) = parser.parse_args()
   if not path.exists(options.outDir):
   	makedirs(options.outDir)
   main(options, args)