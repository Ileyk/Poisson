import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import sys
import h5py
import imageio
import glob
import os
import subprocess

def movie(start,end,step):

    output_folder='../scratch/output/'

    # output location
    path_to_frames='./frames/'
    # If folder does not exist, make it
    if not os.path.exists(path_to_frames):
        os.makedirs(path_to_frames)
    # else, remove all its content
    else:
        files = glob.glob(path_to_frames+'*')
        for f in files:
            os.remove(f)

    start=sys.argv[1]
    end=sys.argv[2]
    step=sys.argv[3]

    input_file=output_folder+'w_000000.dat'
    proc = subprocess.Popen(['wc',input_file], stdout=subprocess.PIPE)
    tmp = proc.communicate()
    N1=int(tmp[0].split()[1]) # number of cells in direction 1
    w=np.zeros(N1)
    f=np.zeros(N1)
    fth=np.zeros(N1)
    x=np.zeros(N1)

    # xmin=-0.5
    # xmax= 0.5
    # L=xmax-xmin
    # for i in range(N1):
    #     x[i]=xmin+(xmax-xmin)*(float(i)+0.5)/float(N1)

    # Read x
    input_file=output_folder+'x_'+str(0).zfill(6)+'.dat'
    file=open(input_file)
    for i in range(N1):
        tmp = file.readline()
        x[i]=float(tmp.split()[0])
    file.close()
    xmin=np.min(x)
    xmax=np.max(x)
    L=xmax-xmin

    # Read rho
    input_file=output_folder+'w_'+str(0).zfill(6)+'.dat'
    file=open(input_file)
    for i in range(N1):
        tmp = file.readline()
        w[i]=float(tmp.split()[0])
    file.close()

    # Read f analytic
    input_file=output_folder+'fth_'+str(0).zfill(6)+'.dat'
    file=open(input_file)
    for i in range(N1):
        tmp = file.readline()
        fth[i]=float(tmp.split()[0])
    file.close()

    for it in range(int(start),int(end)+1,int(step)):

        print("=====",int(100.*((float(it)-float(start)+1)/(float(end)-float(start)+1))),"%","=====")

        input_file=output_folder+'f_'+str(it).zfill(6)+'.dat'
        file=open(input_file)
        for i in range(N1):
            tmp = file.readline()
            f[i]=float(tmp.split()[0])
        file.close()

        fig, axs = plt.subplots(ncols=1,nrows=2)

        ax=axs[0]
        plt.subplots_adjust(wspace=0.,hspace=0.2)
        ax.plot(x,f,'o-',lw=1,color='k',alpha=0.7)
        ax.plot(x,fth,'--',lw=1,color='k')
        # ax.plot(x,-1./np.abs(x),'--',lw=1,color='k')
        plt.text(0.05,0.928,str(it),horizontalalignment='left',
             verticalalignment='center', transform = ax.transAxes, fontsize=14, color='k',bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2'))
        ax.set_xlim(xmin-0.05*L,xmax+0.05*L)
        # ax.set_ylim(-1.1*np.max(np.abs(fth)),1.1*np.max(np.abs(fth)))
        ax.set_ylabel('Function (solid) and analytical (dashed)')
        # ax.set_xlim(0.,27.)

        ax=axs[1]
        rel=(f-fth)/fth # relative difference
        plt.subplots_adjust(wspace=0.,hspace=0.2)
        ax.plot(x,rel,'o-',lw=1,color='k',alpha=0.7)
        # ax.plot(x,np.zeros(len(x)),'-',lw=1,color='k')
        ax.hlines(0.,xmin-0.1*L,xmax+0.1*L,linestyles ="-",colors ='k',lw=1)
        # ax.set_ylim(-1.1*np.max(np.abs(rel)),1.1*np.max(np.abs(rel)))
        ax.set_xlim(xmin-0.05*L,xmax+0.05*L)
        ax.set_ylim(-0.0001,0.0001)
        ax.set_ylabel('Relative difference')

        fig.set_size_inches(18,9)
        plt.savefig(path_to_frames+'frame_'+str(it).zfill(8)+'.tif',format='tif',dpi=150,bbox_inches='tight')
        plt.close()

    filename_gif  =path_to_frames+'movie.gif'
    with imageio.get_writer(filename_gif, mode='I',duration=0.1) as writer:
        filelist = glob.glob(os.path.join(path_to_frames, '*.tif'))
        for filename in sorted(filelist):
            image = imageio.imread(filename)
            writer.append_data(image)

movie(sys.argv[1],sys.argv[2],sys.argv[3])
