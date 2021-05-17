# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 20:15:43 2021

@author: aserafeim
"""

import pandas as pd
from scipy.fftpack import fft, fftfreq, ifft
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

N_sheets=9
plt.close('all')
elements=['Mn Wt%','Si Wt%']

class lineFFT:
    
    def __init__(self,filename,N_sheets,stepsize=5,elements=['Mn Wt%','Si Wt%'],threshold=0.1):
        self.filename=filename
        self.N_sheets=N_sheets
        self.L=stepsize
        self.elements=elements
        self.threshold=threshold
        
        
    def open_clean(self):
        df = pd.read_excel (self.filename,[i for i in range(N_sheets)])
        
        for i in range(len(df)):           
            ind=np.where(df[i][self.elements[0]]>7.5)
            df[i][self.elements[0]+'_C']=df[i][self.elements[0]]
            df[i][self.elements[0]+'_C'][ind[0]]=df[i][self.elements[0]].mean()
            
        return df
#        df_norm = pd.read_exciel (r'Lines_10pct.xlsx',[i for i in range(N_sheets)])
        
    
    def calc_fft(self):
        
        d=self.open_clean() # df contains the composition after inverse FFT
        
        Lf_real={} #Dictionary with the segregation lengths in real space
        ft_mn_v={} #Dictionary with the values after FFT
        Lf_v={} # Dictionary with segregation length in frequency space
        thres_2=0.5
        x_new, x_temp=[], []
        Lf_conc_temp=[]
        Lf_conc={} # Dictionary with segreagation length in real space after all of the simmilar frequencies have been averaged
        power=[]
        for i in range(self.N_sheets):
            N=len(df[i])
            #FFT
            ft_mn=fft(np.array(df[i][self.elements[0]+'_C']))
            ft_mn_v[i]=ft_mn[1:N//2] 
            #Calculate frequencies
            Lf=fftfreq(N, self.L)
            Lf_v[i]=Lf[1:N//2]
            #Filter out frequencies based on the 5% highest
            sort_ft=np.sort((2/N)*np.abs(ft_mn))
            #Make an index vector to find which values of ft are higher than the cutoff value of 10%
            index=(2/N)*np.abs(ft_mn)>sort_ft[-round(self.threshold*len(sort_ft))]
            #Zero out the low frequencies
            ft_clean=ft_mn*index
            power.append(sum(abs(ft_clean))/sum(abs(ft_mn)))
            #Calculate inverse frequncies
            Lf_clean=(1/Lf)*index
            Lf_clean=Lf_clean[1:N//2] #take half the elements due to symmetry
            #Assigng values to a dictoíonary
            Lf_real[i]=Lf_clean
            #Inverse FFT
            y=ifft(ft_clean)
            #Addting the results to the pandas datastructure
            df[i]['FFT '+elements[0]]=y
            
            #Average similar segregation length based on thres_2
            x_new, x_temp=[], []
            Lf_conc_temp=[]
            Lf_nozeros=Lf_real[i][np.nonzero(Lf_real[i])]
            x_temp.append(Lf_nozeros[0])
            print ('line'+str(i))

            
            for j in range(1,len(Lf_nozeros),1):
                # print((Lf_nozeros[j]-x_temp[0])/Lf_nozeros[j])
                if (np.abs(Lf_nozeros[j]-x_temp[0]))/Lf_nozeros[j]<thres_2:
                    x_temp.append(Lf_nozeros[j])
                else:
                    Lf_conc_temp.append(stat.mean(x_temp))
                    x_temp=[]
                    x_temp.append(Lf_nozeros[j])
                
                if j==len(Lf_nozeros)-1:
                    Lf_conc_temp.append(Lf_nozeros[j])
            print(Lf_conc_temp)
            Lf_conc[i]=Lf_conc_temp
            
        
        return df, Lf_real, ft_mn_v, Lf_v, Lf_conc, power
    # def average_per_line():
    #     for i in range






def subplotting(df,Lf_real,elements,name):
    i=0
    for i in range(len(df)):
        fig, ax = plt.subplots(3,1) 
        ax[0].plot(df[i]['Distance '], df[i]['FFT '+elements[0]])
        ax[0].plot(df[i]['Distance '],df[i]['Mn Wt%'])
        ax[0].set_ylabel('Mn Concentration')
        ax[0].set_xlabel('Distance (μm)')
        
        
        ax[1].plot(df[i]['Distance '], df[i]['FFT '+elements[0]])
        ax[1].set_ylabel('Mn Concentration')
        ax[1].set_xlabel('Distance (μm)')
        
        ax[2].plot(Lf_real[i])
        ax[2].set_yscale('log')
        ax[2].set_ylabel('Band separation (μm)')
        ax[2].set_xlabel('Harmonic number')
        fig.savefig('Subplot_'+str(i)+'_'+elements[0]+'.png',dpi=600,bbox_inches='tight')

def plotfreq(ft,Lf,thres):
    #Plotting frequency-magnitude diagrams for all lines
    for i in range(len(Lf)):
        fig, ax =plt.subplots()
        # plt.figure()
        ax.plot(Lf[i][1:],np.abs(ft[i][1:]))
        ax.set_xlabel('Frequency (1/μm)',fontsize=16)
        ax.set_ylabel('Signal magnitude',fontsize=16)     
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(16)
        fig.savefig('MagnitudevsFrequency'+'_'+elements[0]+'Thresh'+str(thres)+'line'+str(i)+'.png',dpi=600,bbox_inches='tight')

def adjust_lengths(Lf_conc,N_freq):
    #dictionary with each component and a list of the freq of each line for this component.
    comp_L_conc={}
    final_L={}
    for j in range(N_freq):
        comp_L_conc[j]=[]
        final_L[j]=[]
    for j in range(N_freq):
        for i in range(len(Lf_conc)):
                comp_L_conc[j].append(Lf_conc[i][j])
    for j in range(N_freq):
        final_L[j].append(np.mean(comp_L_conc[j][0:2]))
        final_L[j].append(np.mean(comp_L_conc[j][3:5]))
        final_L[j].append(np.mean(comp_L_conc[j][6:8]))
    return comp_L_conc, final_L
            
        
        

# def twoline(Lf1,Lf2):
    
file20='Lines_20pct.xlsx'
file10='Lines_10pct.xlsx'
thres=0.5
line20=lineFFT(file20,9,threshold=thres)

df=line20.open_clean()


line10=lineFFT(file10,8,threshold=thres)


y , Lf_real,ft_mn,Lf, Lf_conc, power =line20.calc_fft()

comp_L_conc,final_L=adjust_lengths(Lf_conc,5)
# y50 , Lf_real50,ft_mn50,Lf50, Lf_conc50, power50 =line10.calc_fft()
# subplotting(y,Lf_real,elements,file20)
# 

plotfreq(ft_mn,Lf,thres)

# fig, ax =plt.subplots()
# ax.plot(Lf_conc[4], label='10 pct threshold')
# ax.plot(Lf_conc50[4],label='50 pct threshold')
# ax.legend()

# fig2, ax =plt.subplots()
# ax.plot(power, label='10 pct threshold')
# ax.plot(power50,label='50 pct threshold')
# ax.legend()

# fig3, ax =plt.subplots()
# ax.plot(Lf_real[4], label='10 pct threshold')
# ax.plot(Lf_real50[4],label='50 pct threshold')
# ax.legend()
    
# fig,ax=plt.subplots()
# for i in range(3):
#     ax.plot(Lf_real[i])

fig2, ax =plt.subplots()

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
 	label.set_fontsize(12)
    
sections=['Top', 'Middle', 'Center']
ax.plot(sections,final_L[2],label='Component 2')
ax.plot(sections,final_L[3],label='Component 3')
ax.plot(sections,final_L[4],label='Component 4')
ax.set_xlabel('Section')
ax.set_ylabel('Segregation lenght (μm)')
ax.legend()
fig2.savefig('segregation_lengths'+'_'+elements[0]+'Thresh'+str(thres)+'.png',dpi=600,bbox_inches='tight')
    
    


