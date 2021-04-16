# vim: set et sw=4 ts=4 nu fdm=indent:
# coding: utf8


import os
import numpy as np
import matplotlib.pyplot as plt
    

def MOKE(system, N=100, theta=0, Bmax=1.0, Ntest=0, plotname="", transverse=False, detailed=False, 
         animated=False, anim_name="kerr", anim_size=10, anim_center=[0, 0],
         filename=""):
    
    N = int(N/2.0)
    Bs1 = (np.linspace(0, 2*N, N,endpoint=False)-N)*Bmax/N
    Bs2 = (N-np.linspace(0, 2*N,N+1,endpoint=True))*Bmax/N
    Bs = np.concatenate((Bs1, Bs2))
    
    if Bmax == 0:
        print("Error: null applied field. MOKE measurement can't be done.")
        return 0, 0
    
    if system.dim == "2D":
        Bdir = [np.cos(theta*np.pi/180), np.sin(theta*np.pi/180)]
    elif system.dim == "3D":
        Bdir = [np.cos(theta*np.pi/180), np.sin(theta*np.pi/180), 0]
    system.Bu = Bdir
    
    magl = []
    magt = []
    
    anim_folder=anim_name+"_animation"
    
    if animated:
        if not os.path.isdir(anim_folder):
            os.mkdir(anim_folder)
            print('\n"{}" folder created.\n'.format(anim_folder))
        print('\nSaving the Kerr cycles frames in the "{}" folder.'.format(anim_folder))
    
    print("\nDoing MOKE measurement for \u03B8={}°... (N={})\n".format(theta, len(Bs)))
    i=0
    for Bapp in Bs:
        system.B = Bapp
        
        
        if Ntest != 0:
            system.find_min(N=Ntest, detailed=True, display=False)
        else:
            system.relax(display=detailed)
        
        if animated:
            pname = system.draw(name="{}_{}".format(anim_name, i), size=anim_size, center=anim_center, display=False)
            os.rename(pname, os.path.join(anim_folder, pname))
            print("Frame {} saved.".format(i))
            i += 1
        
        magl.append(np.dot(system.mag, Bdir)*system.M)
        
        if system.dim == "2D":
            Htdir = [Bdir[1], -Bdir[0]]
        elif system.dim == "3D":
            Htdir = [Bdir[1], -Bdir[0], 0]
        magt.append(np.dot(system.mag, Htdir)*system.M)
    
    print("\nMOKE measurement done.\n")
    
    magl = np.array(magl)/system.Ms
    magt = np.array(magt)/system.Ms
    
    if plotname != "":
        plt.plot(Bs, magl, label="M$_L$")
        
        if transverse:
            plt.plot(Bs, magt, label="M$_T$")
            plt.axhline(y=0, color="black", linestyle="dotted", alpha=0.25)
            
            plt.legend(loc="upper left", fontsize="small", bbox_to_anchor=[0, 0.95])
        
        plt.axvline(x=0, color="black", linestyle="dotted", alpha=0.25)
        
        plt.axhline(y=+1, color="blue", linestyle="dashed", alpha=0.5)
        plt.axhline(y=-1, color="blue", linestyle="dashed", alpha=0.5)
        
        plt.xlabel("Applied field B (T)")
        plt.ylabel("M/M$_s$")
        plt.title("Kerr cycle for \u03B8={}°".format(theta))
        
        plt.savefig(plotname+".png", dpi=200)
        plt.clf()
        print('Kerr cycle saved as "{}.png".'.format(plotname))
    
    '''
    print("")
    rem = 0
    for i in range(len(Hs)):
        if Hs[i]==0:
            rem += abs(magl[i])
            break
    
    if detailed:
        print("\nRemanent magnetization measured for theta={}° :\n{} Ms".format(theta, rem))
    '''
    if filename != "":
        file = open(filename, "w")
        print('\nWriting data in the "{}" file...'.format(filename))
        
        file.write("N\tB\tMl\tMt\n")
        for i in range(len(Bs)):
            file.write("{}\t{}\t{}\t{}\n".format(i, Bs[i], magl[i], magt[i]))
            
        print("Done writing.\n")
        file.close()
        
    return Bs, magl, magt


def MOKE_angular(system, angles, N=100, Hmax=1.0, foldername=""):
    if len(angles) == 0:
        print("\nError : no angles given.\n")
        return 0, 0
    
    every_mags = []
    
    try:
        os.mkdir(foldername)
    except:
        pass
    
    for angle in angles:
        if foldername != "":
            plotname = foldername + "_{}".format(angle)
            filename = foldername + "_{}".format(angle)
        
        Hs, mags, magt = MOKE(system, N=N, theta=angle, Hmax=Hmax, plotname=plotname, filename=filename)
        every_mags.append(mags)
        
        plotname = plotname+".png"
        
        if foldername != "":
            os.rename(plotname, os.path.join(foldername, plotname))
            os.rename(filename, os.path.join(foldername, filename))
    
    return Hs, every_mags
    
    
def MOKE_mean(system, Nsamples, theta=0.0, N=100, Bmax=1.0, foldername=""):
    if Nsamples == 0:
        print("\nError : zero number of samples given.\n")
        return 0, 0

    every_mags = []
    
    try:
        os.mkdir(foldername)
    except:
        pass
    
    for i in range(Nsamples):
        if foldername != "":
            plotname = foldername + "_{}".format(i)
            filename = foldername + "_{}".format(i)
        
        Bs, mags, magt = MOKE(system, N=N, theta=theta, Bmax=Bmax, plotname=plotname, filename=filename)
        system.randomize_angles()
        every_mags.append(mags)
        
        plotname = plotname+".png"
        
        if foldername != "":
            os.rename(plotname, os.path.join(foldername, plotname))
            os.rename(filename, os.path.join(foldername, filename))
    
    
    for i in range(Nsamples):
        plt.plot(Bs, every_mags[i], alpha=0.75)
    
    plt.axvline(x=0, color="black", linestyle="dotted", alpha=0.25)
    plt.axhline(y=+1, color="blue", linestyle="dashed", alpha=0.5)
    plt.axhline(y=-1, color="blue", linestyle="dashed", alpha=0.5)
    
    plt.xlabel("Applied field B (T)")
    plt.ylabel("M/M$_s$")
    plt.title("{} Kerr cycles for \u03B8={}°".format(Nsamples, theta))
    
    allplotname = foldername+"_all.png"
    plt.savefig(allplotname, dpi=200)
    plt.clf()
    print('All {} Kerr cycles saved as "{}".\n'.format(Nsamples, allplotname))
    
    os.rename(allplotname, os.path.join(foldername, allplotname))
    
    
    mean = []
    for j in range(len(mags)):
        y = 0.0
        for i in range(Nsamples):
            y += every_mags[i][j]
        mean.append(y/Nsamples)
    
    plt.plot(Bs, mean)
    plt.axvline(x=0, color="black", linestyle="dotted", alpha=0.25)
    plt.axhline(y=+1, color="blue", linestyle="dashed", alpha=0.5)
    plt.axhline(y=-1, color="blue", linestyle="dashed", alpha=0.5)
    
    plt.xlabel("Applied field B (T)")
    plt.ylabel("M/M$_s$")
    plt.title("Mean Kerr cycle of {} cycles for \u03B8={}°".format(Nsamples, theta))
    
    meanplotname = foldername+"_mean.png"
    plt.savefig(meanplotname, dpi=200)
    plt.clf()
    print('Mean Kerr cycle saved as "{}".\n'.format(meanplotname))
    
    os.rename(meanplotname, os.path.join(foldername, meanplotname))
        
    return Bs, every_mags, mean

















