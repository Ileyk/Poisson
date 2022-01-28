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
from matplotlib.colors import LogNorm, SymLogNorm, Normalize

def rho(it):

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

    N1=42
    N2=42
    N=N1*N2
    w=np.zeros((N2,N1))
    f=np.zeros((N2,N1))
    fth=np.zeros((N2,N1))
    r=np.zeros((N2,N1))
    x=np.zeros(N1)
    y=np.zeros(N2)

    # Read x
    input_file=output_folder+'x_'+str(0).zfill(6)+'.dat'
    file=open(input_file)
    for i in range(N1):
        tmp = file.readline()
        x[i]=float(tmp.split()[0])
    for i in range(N2):
        tmp = file.readline()
        y[i]=float(tmp.split()[0])
    print(x)
    print(y)
    file.close()
    xmin=np.min(x)
    xmax=np.max(x)
    ymin=np.min(y)
    ymax=np.max(y)
    Lx=xmax-xmin
    Ly=ymax-ymin

    X, Y = np.meshgrid(x, y)

    input_file=output_folder+'w_'+str(0).zfill(6)+'.dat'
    file=open(input_file)
    for i in range(N1):
        for j in range(N2):
            tmp = file.readline()
            w[j,i]=float(tmp.split()[0])
    file.close()

    # input_file=output_folder+'f_'+str(it).zfill(6)+'.dat'
    input_file=output_folder+'fth_'+str(0).zfill(6)+'.dat'
    file=open(input_file)
    for i in range(N1):
        for j in range(N2):
            tmp = file.readline()
            fth[j,i]=float(tmp.split()[0])
    file.close()

    for it in range(int(start),int(end)+1,int(step)):

        print("=====",int(100.*((float(it)-float(start)+1)/(float(end)-float(start)+1))),"%","=====")

        input_file=output_folder+'f_'+str(it).zfill(6)+'.dat'
        file=open(input_file)
        for i in range(N1):
            for j in range(N2):
                tmp = file.readline()
                f[j,i]=float(tmp.split()[0])
        file.close()

        input_file=output_folder+'r_'+str(it).zfill(6)+'.dat'
        file=open(input_file)
        for i in range(N1):
            for j in range(N2):
                tmp = file.readline()
                r[j,i]=float(tmp.split()[0])
        file.close()

        fig, axs = plt.subplots(ncols=2,nrows=2)

        ax=axs[0,0]
        cax=ax.pcolormesh(X,Y,w,cmap='gist_rainbow')
        cbar=plt.colorbar(cax,ax=ax)
        ax.set_xlim(xmin-0.05*Lx,xmax+0.05*Lx)
        ax.set_ylim(ymin-0.05*Ly,ymax+0.05*Ly)
        ax.title.set_text('Source')

        ax=axs[0,1]

        # plt.subplots_adjust(wspace=0.,hspace=0.2)
        cax=ax.pcolormesh(X,Y,f,cmap='gist_rainbow')
        cbar=plt.colorbar(cax,ax=ax)
        ax.set_xlim(xmin-0.05*Lx,xmax+0.05*Lx)
        ax.set_ylim(ymin-0.05*Ly,ymax+0.05*Ly)
        # ax.set_ylim(-1.1*np.max(np.abs(fth)),1.1*np.max(np.abs(fth)))
        # ax.set_xlim(0.,27.)
        ax.title.set_text('Solution f at iteration '+str(it))

        ax=axs[1,0]

        # plt.subplots_adjust(wspace=0.,hspace=0.2)
        cax=ax.pcolormesh(X,Y,(f-fth)/fth,norm=SymLogNorm(linthresh=1.E-4,linscale=1.E-4,vmin=-1E0,vmax=1E0),cmap='seismic')
        cbar=plt.colorbar(cax,ax=ax)
        ax.set_xlim(xmin-0.05*Lx,xmax+0.05*Lx)
        ax.set_ylim(ymin-0.05*Ly,ymax+0.05*Ly)
        # ax.set_ylim(-1.1*np.max(np.abs(fth)),1.1*np.max(np.abs(fth)))
        # ax.set_xlim(0.,27.)
        ax.title.set_text('(f-fth)/fth')

        ax=axs[1,1]

        cax=ax.pcolormesh(X,Y,r,norm=SymLogNorm(linthresh=1.E-4,linscale=1.E-4,vmin=-1E0,vmax=1E0),cmap='seismic')
        cbar=plt.colorbar(cax,ax=ax)
        ax.set_xlim(xmin-0.05*Lx,xmax+0.05*Lx)
        ax.set_ylim(ymin-0.05*Ly,ymax+0.05*Ly)
        ax.title.set_text('Numerical residual')

        # ax=axs[1]
        # print(np.max(np.abs(f-fth)))
        # rel=(f-fth)/np.abs(fth) # relative difference
        # plt.subplots_adjust(wspace=0.,hspace=0.2)
        # ax.plot(x,rel,'o-',lw=1,color='k',alpha=0.7)
        # # ax.plot(x,np.zeros(len(x)),'-',lw=1,color='k')
        # ax.hlines(0.,xmin-0.1*L,xmax+0.1*L,linestyles ="-",colors ='k',lw=1)
        # # ax.set_ylim(-1.1*np.max(np.abs(rel)),1.1*np.max(np.abs(rel)))
        # ax.set_xlim(xmin-0.05*L,xmax+0.05*L)
        # ax.set_ylim(-0.01,0.01)
        # ax.set_ylabel('Relative difference')

        # fig.set_size_inches(9,9)
        # plt.savefig('frame_'+str(0).zfill(8)+'.tif',format='tif',dpi=150,bbox_inches='tight')
        # plt.close()

        fig.set_size_inches(14,14)
        plt.savefig(path_to_frames+'frame_'+str(it).zfill(8)+'.tif',format='tif',dpi=150,bbox_inches='tight')
        plt.close()

    filename_gif  =path_to_frames+'movie.gif'
    with imageio.get_writer(filename_gif, mode='I',duration=0.1) as writer:
        filelist = glob.glob(os.path.join(path_to_frames, '*.tif'))
        for filename in sorted(filelist):
            image = imageio.imread(filename)
            writer.append_data(image)

rho(sys.argv[1])
