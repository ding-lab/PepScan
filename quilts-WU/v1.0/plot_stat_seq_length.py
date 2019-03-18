#!C:\Python27\python.exe -u
##!cd C:\Users\fenyo\Desktop\xeno
from pandas import Series, DataFrame
import pandas as pd
import matplotlib
matplotlib.use('agg') 
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import numpy
import sys, os
import pylab as P

listing = os.listdir('.')
for infile in listing:
	if infile.find('fasta-length.stat')>=0:
		print infile
		df = pd.read_table(infile)
		types={}
		for column in df.columns[1:]:
			if column.split('-')[-1] in types:
				types[column.split('-')[-1]]+=1
			else:
				types[column.split('-')[-1]]=1

		for column in df.columns[1:]:
			#print column,df[column].max()
			if df[column].max()<=0:
				print 'Empty: '+column

		for type in types:
			print type
			df_=df
			df_=df_.fillna(0)
			fig, (ax1) = plt.subplots(1,figsize=(6,6))
			count=0
			for column in df.columns[1:]:
				if column.split('-')[-1]==type:
					if df[column].max()>0:
						#print column
						plt.plot(df[column],color=P.cm.hot(count/1.5/types[column.split('-')[-1]]))
						ax1.set_yscale('log')
						ax1.set_xscale('log')
					count+=1
			plt.savefig(infile+'-'+type+'-all.png',dpi=72,bbox_inches='tight')
			
		plt.clf()



