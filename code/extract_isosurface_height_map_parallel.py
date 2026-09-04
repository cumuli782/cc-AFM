import sys
import math as M
import time
import argparse
import os
import numpy as np

from multiprocess import shared_memory, Process, Manager, cpu_count, set_start_method, Pool, Lock, Value

lock=Lock()

no_iso=False
only_iso=False

def dist(X1,X2):
    summ=0
    for i in range(len(X1)):
        summ += (X1[i]-X2[i])**2
    summ = summ**(0.5)
    return summ


def d1_to_d3_arr(arr,nx,ny,nz):   #converts 1d density array in 3d array according to the grid parameters to be processed
    if isinstance(arr[0], float):
        d3_arr = np.ndarray(shape=(nx,ny,nz),dtype=np.single)
        
        for i, it in enumerate(arr):
            iz = i // (nx*ny)
            rst= i % (nx*ny)
            iy = rst // nx
            rst = rst % nx
            ix = rst

            d3_arr[ix,iy,iz] = it

    else:
        d3_arr=np.ndarray(shape=(nx,ny,nz,3),dtype=np.single)
        for i, it in enumerate(arr):
            iz = i // (nx*ny)
            rst= i % (nx*ny)
            iy = rst // nx
            rst = rst % nx
            ix = rst

            d3_arr[ix,iy,iz,0] = it[0]
            d3_arr[ix,iy,iz,1] = it[1]
            d3_arr[ix,iy,iz,2] = it[2]

    return d3_arr



def linearzerlegung(pos,vecs,nx,ny,nz,offset=[0,0,0]): #Transforms position into new coordinate system (AFM grid to chg grid or vice versa)
    mat = [[vecs[0][0],vecs[1][0],vecs[2][0]],[vecs[0][1],vecs[1][1],vecs[2][1]],[vecs[0][2],vecs[1][2],vecs[2][2]]]
    #print(mat)
    #print("Input_pos:",pos)
    wrk_pos = [pos[i] -offset[i] for i in range(3)]
    
    if vecs[0][0] == 0.0 or vecs[1][1] == 0.0:
        tp_mat = [mat[1],mat[0]]
        mat = tp_mat
        tp_pos = [pos[1],pos[0]]
        wrk_pos = tp_pos

    fac = -mat[1][0] / mat[0][0]
    #mat[1][0] += mat[0][0]*fac
    mat[1][1] += mat[0][1]*fac
    wrk_pos[1] += wrk_pos[0]*fac

    #mat[0][2] muss 0 sein, daher nichts weiter noetig
    #mat[2][0] ist auch schon 0

    #print ("fac:",fac)
    #mat[2][1] ist 0, kein Schritt hierfuer erforderlich
    #ebenso fuer mat[1][2] und mat[0][2]
    
    fac = -mat[0][1] / mat[1][1]
    #mat[0][1] += mat[1][1]*fac
    #mat[0][0] += mat[1][0]*fac
    wrk_pos[0] += wrk_pos[1]*fac

    #print ("fac",fac)
    #mat 

    #print(wrk_pos[1])
    wrk_pos[0] /= mat[0][0]
    wrk_pos[1] /= mat[1][1]
    wrk_pos[2] /= mat[2][2]

    #print(wrk_pos[1])
    
    idx_x = wrk_pos[0]*nx
    idx_y = wrk_pos[1]*ny
    idx_z = wrk_pos[2]*nz

    return idx_x, idx_y, idx_z

def interpolate_probe(ix,iy,iz,d3_arr):#Interpolates the deflections of the probe particle 

    no_probe_interpolation= False

    if no_probe_interpolation:
        ix = int(ix+0.5)
        iy = int(iy+0.5)
        iz = int(iz+0.5)
        return d3_arr[ix,iy,iz]
    
    dist_cutoff=0.5

    #finds closest point
    
    ix_c = int(ix+0.5)
    iy_c = int(iy+0.5)
    iz_c = int(iz+0.5)

    #applies boundary conditions
    
    if ix_c >= len(d3_arr):
        ix_c -= len(d3_arr)

    if iy_c >= len(d3_arr[0]):
        iy_c -=	len(d3_arr[0])
    if iz_c >= len(d3_arr[0,0]):
        iz_c -=	len(d3_arr[0,0])

    #finds next and second nearest point in x, y and z direction
    
    if ix_c == int(ix):
        ix_2nd = ix_c + 1
        if ix_c == len(d3_arr) -1:
            ix_2nd = 0
        ix_3rd = ix_c - 1
    
    else:
        ix_2nd = ix_c - 1
        ix_3rd = ix_c + 1
        if ix_c == len(d3_arr) -1:
            ix_3rd = 0
    

    if iy_c == int(iy):
        iy_2nd = iy_c + 1
        if iy_c == len(d3_arr[0]) -1:
            iy_2nd = 0
        iy_3rd = iy_c - 1
    
    else:
        iy_2nd = iy_c - 1
        iy_3rd = iy_c + 1
        if iy_c == len(d3_arr[0]) -1:
            iy_3rd = 0
    

    if iz_c == int(iz):
        iz_2nd = iz_c + 1
        if iz_c == len(d3_arr[0,0]) -1:
            iz_2nd = 0
        iz_3rd = iz_c -1
    
    else:
        iz_2nd = iz_c - 1
        iz_3rd = iz_c +1
        if iz_c == len(d3_arr[0,0])-1:
            iz_3rd = 0
    

    #determines from which side the interpolation takes place. If both next points are to far away, no interpolation is done
            
    xcase=""
    if dist(d3_arr[ix_c,iy_c,iz_c],d3_arr[ix_2nd,iy_c,iz_c]) < dist_cutoff:
        xcase = "c-2nd"
    elif dist(d3_arr[ix_c,iy_c,iz_c],d3_arr[ix_3rd,iy_c,iz_c]) < dist_cutoff:
        xcase = "c-3rd"
    else:
        xcase = "c-c"

    ycase=""
    if dist(d3_arr[ix_c,iy_c,iz_c],d3_arr[ix_c,iy_2nd,iz_c]) < dist_cutoff:
        ycase =	"c-2nd"
    elif dist(d3_arr[ix_c,iy_c,iz_c],d3_arr[ix_c,iy_3rd,iz_c]) < dist_cutoff:
        ycase =	"c-3rd"
    else:
        ycase =	"c-c"

    zcase=""
    if dist(d3_arr[ix_c,iy_c,iz_c],d3_arr[ix_c,iy_c,iz_2nd]) < dist_cutoff:
        zcase =	"c-2nd"
    elif dist(d3_arr[ix_c,iy_c,iz_c],d3_arr[ix_c,iy_c,iz_3rd]) < dist_cutoff:
        zcase =	"c-3rd"
    else:
        zcase =	"c-c"

    # Interpolates the position
    dX= []
    if xcase == "c-2nd":
        di = abs(ix - ix_c)
        dX = [(d3_arr[ix_2nd,iy_c,iz_c,i]- d3_arr[ix_c,iy_c,iz_c,i]) * di for i in range(3) ] 
   
    elif xcase == "c-3rd":
        di = abs(ix-ix_c)
        dX = [(d3_arr[ix_c,iy_c,iz_c,i] - d3_arr[ix_3rd,iy_c,iz_c,i]) * di for i in range(3) ]

    elif xcase == "c-c":
        dX = [0,0,0]

    dY =[]
    if ycase == "c-2nd":
        di = abs(iy - iy_c)
        dY = [(d3_arr[ix_c,iy_2nd,iz_c,i]- d3_arr[ix_c,iy_c,iz_c,i]) * di for i in range(3) ]

    elif ycase == "c-3rd":
        di = abs(iy-iy_c)
        dY = [(d3_arr[ix_c,iy_c,iz_c,i] - d3_arr[ix_c,iy_3rd,iz_c,i]) * di for i in range(3) ]

    elif ycase == "c-c":
        dY = [0,0,0]


    dZ =[]
    if zcase == "c-2nd":
        di = abs(iz - iz_c)
        dZ = [(d3_arr[ix_c,iy_c,iz_2nd,i]- d3_arr[ix_c,iy_c,iz_c,i]) * di for i in range(3) ]

    elif zcase == "c-3rd":
        di = abs(iz-iz_c)
        dZ = [(d3_arr[ix_c,iy_c,iz_c,i] - d3_arr[ix_c,iy_c,iz_3rd,i]) * di for i in range(3) ]

    elif zcase == "c-c":
        dZ = [0,0,0]

    try:
        probe_pos = [ d3_arr[ix_c,iy_c,iz_c,i] + dX[i] + dY[i] + dZ[i] for i in range(3) ]
    
    except IndexError:
        print("Index Error in interpolate probe:",flush=True)
        print("not interpolated deflection:", d3_arr[ix_c][iy_c][iz_c], "dX:", dX, "dY:", dY, "dZ:", dZ)
        return [0,0,0], "Probe interpolation Error"
        
    return probe_pos, 0

def interpolate_charge(ix,iy,iz,d3_arr):#Interpolates the charge measured at the probe particle position
    try:
        dr_x = (ix - int(ix)) * (d3_arr[int(ix+1),int(iy),int(iz)] - d3_arr[int(ix),int(iy),int(iz)])
    except IndexError:#Applies boundary condition if necessary
        flag = True
        if ix >= len(d3_arr) -1:
            ix -= len(d3_arr)
            flag = False


        if iy >= len(d3_arr[0]) -1:
            iy -= len(d3_arr[0])
            flag = False

            

        if iz >= len(d3_arr[0,0]) -1:
            iz -= len(d3_arr[0,0])
            flag = False

            
        if flag:
             print("Error while interpolating charge in x direction, indices", ix, iy, iz)
             return 0, "Charge interpolation Error, possible grid mismatch"
             #sys.exit()

        dr_x = (ix - int(ix)) * (d3_arr[int(ix+1),int(iy),int(iz)] - d3_arr[int(ix),int(iy),int(iz)])
        
    try:
        dr_y = (iy - int(iy)) * (d3_arr[int(ix),int(iy+1),int(iz)] - d3_arr[int(ix),int(iy),int(iz)])
    except IndexError:
        flag=True
        if ix >= len(d3_arr) -1:
            ix -= len(d3_arr)
            flag=False


        
        if iy >= len(d3_arr[0]) -1:
            iy -= len(d3_arr[0])
            flag=False


        if iz >= len(d3_arr[0,0]) -1:
            iz -= len(d3_arr[0,0])
            flag=False

        if flag:   
             print("Error while interpolating charge in y direction, indices", ix, iy, iz)
             return 0, "Charge interpolation Error, possible grid mismatch"
             #sys.exit()

        dr_y = (iy - int(iy)) * (d3_arr[int(ix),int(iy+1),int(iz)] - d3_arr[int(ix),int(iy),int(iz)])
    
    try:         
        dr_z = (iz - int(iz)) * (d3_arr[int(ix),int(iy),int(iz+1)] - d3_arr[int(ix),int(iy),int(iz)])
    except IndexError:
        flag=True
        if ix >= len(d3_arr) -1:
            ix -= len(d3_arr)
            flag=False



        if iy >= len(d3_arr[0]) -1:
            iy -= len(d3_arr[0])
            flag=False

        
        if iz >= len(d3_arr[0,0]) -1:
            iz -= len(d3_arr[0,0])
            flag=False

        if flag:
             print("Error while interpolating charge in z direction, indices", ix, iy, iz)
             return 0, "Charge interpolation Error, possible grid mismatch"
             #sys.exit()

        dr_z = (iz - int(iz)) * (d3_arr[int(ix),int(iy),int(iz+1)] - d3_arr[int(ix),int(iy),int(iz)])

        
    rho = d3_arr[int(ix),int(iy),int(iz)] + dr_x + dr_y + dr_z
    return rho, 0



def pp_t_dist(R,pos_r,pos_ur,x,y,z):#calculates the difference between the distances between the probe-particle and the tip in the current and equilibrium position
    
    dX = [pos_r[i] -pos_ur[i] for i in range(3)]

    if dX[0] > R:
        dX = [dX[i]-x[i] for i in range(3)]
    if dX[0] > R:
        dX = [dX[i]-y[i] for i in range(3)]
    if dX[1] > R:
        dX = [dX[i]-y[i] for i in range(3)]
    if dX[1] > R:
        dX = [dX[i]-x[i] for i in range(3)]
    if dX[0] < -R:
        dX = [dX[i]+x[i] for i in range(3)]
    if dX[0] < -R:
        dX = [dX[i]+y[i] for i in range(3)]
    if dX[1] < -R:
        dX = [dX[i]+y[i] for i in range(3)]
    if dX[1] < -R:
        dX = [dX[i]+x[i] for i in range(3)]
    
    dR =  (dX[0]**2 + dX[1]**2 + (-R + dX[2])**2)** 0.5 - R    
    
    
    return dR

def charge_at_unrelaxed_pp_pos(d3_arr,chg_arr_info,pp_d3_arr,deflection_arr_info, deflection_offset_idx,pp_ix,pp_iy,pp_iz, tip_flipped,R,b_decay,tip_flip_cutoff):#gets the charge density at the position of the relaxed probe particle (function name inconvenient)

    #extracting parameters
    nx = chg_arr_info[1]
    ny = chg_arr_info[2]
    nz = chg_arr_info[3]
    x = chg_arr_info[4]
    y = chg_arr_info[5]
    z = [0,0,chg_arr_info[6]]


    pp_nx = deflection_arr_info[1]
    pp_ny = deflection_arr_info[2]
    pp_nz = deflection_arr_info[3]
    pp_x = deflection_arr_info[4]
    pp_y = deflection_arr_info[5]
    pp_z = [0,0,deflection_arr_info[6]]
    pp_offset = deflection_arr_info[7]


    #Getting unrelaxed probe-particle position

    pp_unrelaxed = [pp_ix/pp_nx * pp_x[i] + pp_iy/pp_ny * pp_y[i] + pp_iz/pp_nz * pp_z[i] + pp_offset[i] for i in range(3)]

    #Getting relaxed probe-particle position
    pp_pos, errorcode = interpolate_probe(pp_ix,pp_iy,pp_iz + deflection_offset_idx,pp_d3_arr)
    if errorcode != 0:
        return 0, 0, errorcode
    pp_pos[2] -= deflection_offset_idx * pp_z[2] / pp_nz

    if pp_unrelaxed[2] + R <= pp_pos[2] and tip_flip_cutoff: #Check if the probe-particle was bend higher than the tip, if yes calculate density at tip instead
        new_pos = pp_unrelaxed.copy()
        new_pos[2] += R
        dens_ix, dens_iy, dens_iz = linearzerlegung(new_pos,[x,y,z],nx,ny,nz)
        chg_val, errorcode = interpolate_charge(dens_ix, dens_iy, dens_iz, d3_arr)
        if errorcode !=0:
            return 0, 0, errorcode
        tip_flipped +=1

        
        return chg_val, tip_flipped, 0
    

    #Calculate probe-particle -tip distance and the resulting decay in tunneling current
    probe_tip_ddist= pp_t_dist(R,pp_pos,pp_unrelaxed,x,y,z)
    try:
        probe_tip_decay = M.exp(-b_decay * probe_tip_ddist)
    except OverflowError:
        print("Overflow Error. probe_tip_ddist =",probe_tip_ddist)
        #chg_arr_mem.close()
        #pp_arr_mem.close()
        return 0, 0, "Error pp_decay overflow"
        #sys.exit()

    #calculate the density at the relaxed probe-particle position
    dens_ix, dens_iy, dens_iz = linearzerlegung(pp_pos,[x,y,z],nx,ny,nz)
    chg_val, errorcode = interpolate_charge(dens_ix, dens_iy, dens_iz, d3_arr)
    if errorcode != 0:
        return 0, 0, errorcode

    #sclae charge with tunneling decay
    chg_val *= probe_tip_decay


    return chg_val, tip_flipped, 0


def average_density(dens_amp_arr):#Average charges over the set amplitude
    chg=0
    for it in dens_amp_arr:
        chg += it
    chg /= len(dens_amp_arr)
    return chg

def check_iso_at_deflection(ix,iy,iz,chg_arr_info,deflection_arr_info,thresh,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp):#Scan z column for the charge iso value, including deflections and averaging. Starts at z value provided by scan_z

    height_interpolation = True

    #Extracts parameters
    nx = chg_arr_info[1]
    ny = chg_arr_info[2]
    nz = chg_arr_info[3]
    x = chg_arr_info[4]
    y = chg_arr_info[5]
    z = [0,0,chg_arr_info[6]]

    #load shared memory objects
    chg_arr_mem= shared_memory.SharedMemory(name=chg_arr_info[0])
    d3_arr = np.ndarray((nx,ny,nz),dtype=np.single,buffer=chg_arr_mem.buf)


    pp_nx = deflection_arr_info[1]
    pp_ny = deflection_arr_info[2]
    pp_nz = deflection_arr_info[3]
    pp_x = deflection_arr_info[4]
    pp_y = deflection_arr_info[5]
    pp_z = [0,0,deflection_arr_info[6]]
    pp_offset = deflection_arr_info[7]

    pp_arr_mem = shared_memory.SharedMemory(name=deflection_arr_info[0])
    pp_d3_arr = np.ndarray((pp_nx,pp_ny,pp_nz,3),dtype=np.single,buffer=pp_arr_mem.buf)
    

    n_Amp = int(Amp / pp_z[2] * pp_nz)
    

    tip_flipped=0
    flipped_flag = False


    #get position from indices
    pos = [ix/nx*x[i] + iy/ny*y[i] + iz/nz * z[i] for i in range(3)]

    
    pp_ix, pp_iy, pp_iz = linearzerlegung(pos,[pp_x,pp_y,pp_z],pp_nx,pp_ny,pp_nz,pp_offset)    
    

    unrelaxed_pos= [pp_ix/pp_nx * pp_x[i] + pp_iy/pp_ny * pp_y[i] + pp_iz/pp_nz * pp_z[i] + pp_offset[i] for i in range(3)]
    
    unrelaxed_pos_save = unrelaxed_pos.copy()
    

    
    #if pp_ix < 0 or pp_iy < 0:
    #    print("Das kann nicht funktionieren, pp_ix, pp_iy2 = ", pp_ix, pp_iy)

    #get charges inside the oscillation amplitude
    dens_amp_array=[]

    
    for i in range(-n_Amp, n_Amp + 1):
        if pp_ix < 0 or pp_iy < 0 or pp_ix >= pp_nx or pp_iy >= pp_ny:
            pos2= pos
            pos2[2] += i * pp_z[2]/pp_nz
            ix2, iy2, iz2 = linearzerlegung(pos2,[x,y,z],nx,ny,nz)
            chg, errorcode = interpolate_charge(ix2,iy2,iz2,d3_arr)
            tip_flipped = 0
        else:    
            chg, tip_flipped, errorcode = charge_at_unrelaxed_pp_pos(d3_arr,chg_arr_info,pp_d3_arr,deflection_arr_info, deflection_offset_idx,pp_ix,pp_iy,pp_iz + i, tip_flipped,R,b_decay,tip_flip_cutoff)
        if errorcode != 0:#Close shared memory acess after error
            chg_arr_mem.close()
            pp_arr_mem.close()
            return 0, [0,0,0], [0,0,0], errorcode
        if tip_flipped >0:
            chg = 100
            flipped_flag ==True
            tip_flipped = 0
        dens_amp_array.append(chg)
        
    
    new_chg_val = average_density(dens_amp_array)#average charges

    
    
    
    if new_chg_val > thresh: #upwards scan for iso value
        chg_old_val=0
        while new_chg_val > thresh:
            # During the upwards scan the probe particle should not flip, or rather it is not relevant, therefore it is not checked here
            if pp_iz + deflection_offset_idx - n_Amp  < 0 or pp_iz + deflection_offset_idx + n_Amp >= pp_nz:
                print("Warning: End of the tip deflections reached (pp_iz =",pp_iz, "), upscan")
                chg_arr_mem.close()
                pp_arr_mem.close()
                return 0, [0,0,0], [0,0,0], "Tip deflection end Error. Check scan range of constant height AFMs"
            chg_old_val = new_chg_val
            pp_iz += 1

            new_chg_val -= dens_amp_array[0]/len(dens_amp_array) #remove charge corresponding to lowest height from average
            for i in range(len(dens_amp_array)-1):
                dens_amp_array[i] = dens_amp_array[i+1]

            #calculate next charge
            if pp_ix < 0 or pp_iy < 0 or pp_ix >= pp_nx or pp_iy >= pp_ny:
                pos2= pos
                pos2[2] =  (pp_iz+n_Amp)/pp_nz * pp_z[2] + pp_offset[2]
                ix2, iy2, iz2 = linearzerlegung(pos2,[x,y,z],nx,ny,nz)
                chg_val, errorcode = interpolate_charge(ix2,iy2,iz2,d3_arr)
                tip_flipped = 0
            else:
                chg_val, tip_flipped, errorcode =charge_at_unrelaxed_pp_pos(d3_arr,chg_arr_info,pp_d3_arr,deflection_arr_info, deflection_offset_idx,pp_ix,pp_iy,pp_iz + n_Amp, tip_flipped,R,b_decay,tip_flip_cutoff)

            if errorcode != 0:
                chg_arr_mem.close()
                pp_arr_mem.close()
                return 0, [0,0,0], [0,0,0], errorcode
            if tip_flipped > 0:
                chg_val = 100
                #print("Achtung: Tip flipped in upward scan!")
                tip_flipped = 0
                
            #include new charge in average
            dens_amp_array[-1] = chg_val
            new_chg_val += dens_amp_array[-1]/len(dens_amp_array)
            
        if height_interpolation:#Interpolate height between the final steps
            diff = new_chg_val - chg_old_val
            dz = (thresh - chg_old_val)/diff
            pp_iz += dz

        new_z = pp_iz / pp_nz * pp_z[2] + pp_offset[2]

    else:
        chg_old_val = 0
        while new_chg_val < thresh: #downward scan for iso value
                    
            if pp_iz + deflection_offset_idx - n_Amp < 0 or pp_iz + deflection_offset_idx + n_Amp >= pp_nz:
                print("Warning: End of tip deflections reached (pp_iz =",pp_iz, "delta dist tip =",probe_tip_ddist, ") downscan")
                chg_arr_mem.close()
                pp_arr_mem.close()
                return 0, [0,0,0], [0,0,0], "Tip deflection end Error. Check scan range of constant height AFMs"

            chg_old_val = new_chg_val
            pp_iz -= 1

            new_chg_val -= dens_amp_array[-1]/len(dens_amp_array)
            for i in range(len(dens_amp_array)-1,0,-1):
                dens_amp_array[i] = dens_amp_array[i-1]

            if pp_ix < 0 or pp_iy < 0 or pp_ix >= pp_nx or pp_iy >= pp_ny:
                pos2= pos
                pos2[2] = (pp_iz-n_Amp)/pp_nz * pp_z[2] + pp_offset[2]
                ix2, iy2, iz2 = linearzerlegung(pos2,[x,y,z],nx,ny,nz)
                chg_val, errorcode = interpolate_charge(ix2,iy2,iz2,d3_arr)
                tip_flipped = 0
            else:
                chg_val, tip_flipped, errorcode = charge_at_unrelaxed_pp_pos(d3_arr,chg_arr_info,pp_d3_arr,deflection_arr_info, deflection_offset_idx,pp_ix,pp_iy,pp_iz - n_Amp, tip_flipped,R,b_decay,tip_flip_cutoff)

            if errorcode !=0:
                chg_arr_mem.close()
                pp_arr_mem.close()
                return 0, [0,0,0], [0,0,0], errorcode
            if tip_flipped >0:
                tip_flipped = 0

            dens_amp_array[0] = chg_val
            
            new_chg_val += dens_amp_array[0]/len(dens_amp_array)

        
            
        if height_interpolation:
            diff = new_chg_val - chg_old_val
            dz = (thresh - chg_old_val)/diff
            pp_iz -= dz

        new_z = pp_iz / pp_nz * pp_z[2] + pp_offset[2]

    
    #Get position of the final deflection
    unrelaxed_pos= [pp_ix/pp_nx * pp_x[i] + pp_iy/pp_ny * pp_y[i] + pp_iz/pp_nz * pp_z[i] + pp_offset[i] for i in range(3)]
    ini_deflection, errorcode= interpolate_probe(pp_ix,pp_iy,pp_iz + deflection_offset_idx,pp_d3_arr)
    if errorcode !=0:
        chg_arr_mem.close()
        pp_arr_mem.close()
        return 0, [0,0,0], [0,0,0], errorcode
    ini_deflection[2] -= deflection_offset_idx * pp_z[2] / pp_nz

    #shift back the height according to the offset
    new_z += deflection_offset_idx /pp_nz * pp_z[2] 
    
    chg_arr_mem.close()
    pp_arr_mem.close()
    return new_z, ini_deflection, unrelaxed_pos, 0
    
    
    
    
    

def scan_z(thresh, z_ax, sstart,ix,iy,chg_arr_info,deflection_arr_info, deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp,shm_height_arr,shm_ini_deflection_array,shm_unrelaxed_pos_array, shm_error_detected):  #scans one z column (col) for the index, where the density is at the provided isovalue (thresh), and returns the corresponding height. Scan is without averaging and deflections, the routine for this is called at the end of scan_z
        
    if (ix*chg_arr_info[2] + iy) % int(chg_arr_info[1]*chg_arr_info[2]/20) ==0:                                                                                                                                                                                                                                                               
        print("%.1f %%" % ((ix*chg_arr_info[2] + iy) / (chg_arr_info[1]*chg_arr_info[2]/100)),flush=True)

    #check if previous run had an error, end run if this is the case
    error_detected_mem = shared_memory.SharedMemory(name=shm_error_detected)
    error_detected = error_detected_mem.buf[0]
        
    if error_detected != 0:
        error_detected_mem.close()
        return 0

    #access shared memory objects
    chg_arr_mem= shared_memory.SharedMemory(name=chg_arr_info[0])
    chg_arr = np.ndarray((chg_arr_info[1],chg_arr_info[2],chg_arr_info[3]), dtype=np.single, buffer=chg_arr_mem.buf)
    col = chg_arr[ix][iy]

    height_arr_mem = shared_memory.SharedMemory(name=shm_height_arr)
    height_arr = np.ndarray((chg_arr_info[1],chg_arr_info[2]), dtype=np.single, buffer=height_arr_mem.buf)

    ini_defl_arr_mem = shared_memory.SharedMemory(name=shm_ini_deflection_array)
    ini_defl_arr = np.ndarray((chg_arr_info[1]*chg_arr_info[2],3), dtype=np.single, buffer=ini_defl_arr_mem.buf)

    unrelaxed_pos_arr_mem = shared_memory.SharedMemory(name=shm_unrelaxed_pos_array)
    unrelaxed_pos_arr = np.ndarray((chg_arr_info[1]*chg_arr_info[2],3), dtype=np.single, buffer=unrelaxed_pos_arr_mem.buf)
    
    nz = len(col)

    istart = M.floor(nz*sstart)   #the scanning starts from this index
    #print(col[istart])
    downsearch = col[istart] < thresh   # Are we already below the threshold? If no, scan downwards
    id_found= "NaN"
    if downsearch:
        id_found = 3*nz
        for i in range(istart, istart-nz,-1):
            if col[i] >= thresh:    # Found the index!
                id_found = i
                drho = -col[i] + col[i+1]   #Interpolation
                di = (thresh - col[i])/drho
                id_found += di
                break
        if id_found < 0:   #Shift the index into the unit cell if necessary
            id_found += nz

    else:# Scan updwards
        id_found = -3*nz
        for i in range(istart, istart + nz):
            iw=i
            if i >= nz:    #Shift index into the unit cell if necessary
                iw -= nz
            if col[iw] < thresh: # Found the index!
                id_found = i
                drho = -col[iw-1] +col[i]    #Interpolation
                di = (thresh - col[i])/drho
                id_found += di
                break
        if id_found > nz:
            id_found -= nz

    height = ""
    if id_found == 3*nz:  #No index found during downwards search: No isosurface in this column!
        height = -100
        no_iso=True
    elif id_found == -3*nz: # No index found during upwards search: column completely within iso surface!
        height = 100
        only_iso=True
    else:
        height = id_found/nz * z_ax #get height from index
        if len(deflection_arr_info) != 0:# If deflections are supplied: Start scan with deflections starting from previous height
            height, ini_deflection, unrelaxed_pos, errorcode = check_iso_at_deflection(ix,iy,id_found,chg_arr_info,deflection_arr_info,thresh,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)
            if errorcode != 0:
                error_detected_mem.buf[0] = 1

            #assign height to shared memory object and close access
            error_detected_mem.close()
            chg_arr_mem.close()
            height_arr[ix][iy] = height
            ini_defl_arr[ix*chg_arr_info[2] + iy][0] = ini_deflection[0]
            ini_defl_arr[ix*chg_arr_info[2] + iy][1] = ini_deflection[1]
            ini_defl_arr[ix*chg_arr_info[2] + iy][2] = ini_deflection[2]
            
            
            unrelaxed_pos_arr[ix*chg_arr_info[2] + iy][0] = unrelaxed_pos[0]
            unrelaxed_pos_arr[ix*chg_arr_info[2] + iy][1] = unrelaxed_pos[1]
            unrelaxed_pos_arr[ix*chg_arr_info[2] + iy][2] = unrelaxed_pos[2]
            
            height_arr_mem.close()
            ini_defl_arr_mem.close()
            unrelaxed_pos_arr_mem.close()
            
            
            return  errorcode

    #assign height to shared memory object and close access (if no tip deflections are given)
    error_detected_mem.close()
    chg_arr_mem.close()

    height_arr[ix][iy] = height
    
    height_arr_mem.close()
    ini_defl_arr_mem.close()
    unrelaxed_pos_arr_mem.close()
    
    return 0

def compute_average_charge(chg_arr_info,deflection_arr_info,deflection_offset_idx, height_arr, outfile):# calculates average charge, only used for testing
    Amp= 1.0 /2

    nx = chg_arr_info[1]
    ny = chg_arr_info[2]
    nz = chg_arr_info[3]

    x= chg_arr_info[4]
    y= chg_arr_info[5]
    
    d3_arr = chg_arr_info[0]
    z = [0,0,chg_arr_info[6]]

    pp_d3_arr = deflection_arr_info[0]
    pp_nx = deflection_arr_info[1]
    pp_ny = deflection_arr_info[2]
    pp_nz = deflection_arr_info[3]
    pp_x = deflection_arr_info[4]
    pp_y = deflection_arr_info[5]
    pp_z = [0,0,deflection_arr_info[6]]
    pp_offset = deflection_arr_info[7]

    outarray=[]
    
    for ix in range(nx):
        outarray.append([])
        for iy in range(ny):
            #unrelaxed PP Position
            pos = [ix/nx * x[i] + iy/ny * y[i] for i in range(3)]
            pos[2]= height_arr[ix][iy] - deflection_offset_idx /pp_nz * pp_z[2]

            #Relaxed PP Position
            pp_ix, pp_iy, pp_iz = linearzerlegung(pos,[pp_x,pp_y,pp_z],pp_nx,pp_ny,pp_nz,pp_offset)
            pp_pos = interpolate_probe(pp_ix,pp_iy,pp_iz + deflection_offset_idx,pp_d3_arr)
            pp_pos[2] -= deflection_offset_idx * pp_z[2]/pp_nz

            #Chg density at relaxed PP
            dens_ix, dens_iy, dens_iz = linearzerlegung(pp_pos,[x,y,z],nx,ny,nz)
            chg= interpolate_charge(dens_ix,dens_iy,dens_iz, d3_arr)

            #Mitteln ueber n_dz Abschnitte innerhalb der Amplitude
            n_dz = 100
            di_z = Amp /pp_z[2] * pp_nz /n_dz

            for i in range(1,n_dz):
                pp_pos_up = interpolate_probe(pp_ix,pp_iy,pp_iz + i * di_z + deflection_offset_idx,pp_d3_arr)
                pp_pos_up[2] -= deflection_offset_idx * pp_z[2]/pp_nz
                dens_ix, dens_iy, dens_iz = linearzerlegung(pp_pos_up,[x,y,z],nx,ny,nz)
                chg += interpolate_charge(dens_ix,dens_iy,dens_iz, d3_arr)

                pp_pos_down = interpolate_probe(pp_ix,pp_iy,pp_iz - i * di_z + deflection_offset_idx,pp_d3_arr)
                pp_pos_down[2] -= deflection_offset_idx * pp_z[2]/pp_nz
                dens_ix, dens_iy, dens_iz = linearzerlegung(pp_pos_down,[x,y,z],nx,ny,nz)
                chg += interpolate_charge(dens_ix,dens_iy,dens_iz, d3_arr)

            chg /= n_dz * 2 - 1
            outarray[-1].append([pos[0],pos[1],pos[2],chg])

    f_out= open(outfile+"_averaged_densities.txt","w")
    for ix in range(nx):
        for iy in range(ny):
            outline =""
            for i in range(4):
                outline += "%f " % (outarray[ix][iy][i])
            outline += "\n"
            f_out.write(outline)

    f_out.close()

            
def compute_height_arr(chg_arr_info, thresh, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp): # Scans the xy-Grid for the provided iso value (thresh) to generate height array
    
    #Get parameters
    nx = chg_arr_info[1]
    ny = chg_arr_info[2]

    z_ax = chg_arr_info[6]
    
    #initialize arrays and corresponding shared memory objects
    
    height_arr = np.zeros(shape=(nx,ny),dtype=np.single)
    ini_deflection_array = np.zeros(shape=(nx*ny,3),dtype=np.single)
    unrelaxed_pos_array = np.zeros(shape=(nx*ny,3),dtype=np.single)

    shm_height_arr= shared_memory.SharedMemory(create=True,size=height_arr.nbytes)
    height_arr_shm=np.ndarray(height_arr.shape,dtype=np.single,buffer=shm_height_arr.buf)
    height_arr_shm[:] = height_arr[:]

    shm_ini_deflection_array= shared_memory.SharedMemory(create=True,size=ini_deflection_array.nbytes)
    ini_deflection_array_shm=np.ndarray(ini_deflection_array.shape,dtype=np.single,buffer=shm_ini_deflection_array.buf)
    ini_deflection_array_shm[:] = ini_deflection_array[:]

    shm_unrelaxed_pos_array= shared_memory.SharedMemory(create=True,size=unrelaxed_pos_array.nbytes)
    unrelaxed_pos_array_shm=np.ndarray(unrelaxed_pos_array.shape,dtype=np.single,buffer=shm_unrelaxed_pos_array.buf)
    unrelaxed_pos_array_shm[:] = unrelaxed_pos_array[:]


    #setup for parallel execution
    print("starting scan",flush=True)
    cpus= cpu_count()

    print("Using", cpus, "Processes")


    shm_error_detected = shared_memory.SharedMemory(create=True,size=1)
    error_detected = shm_error_detected.buf
    error_detected[0] = 0


    workerpool = Pool(cpus)
    errorcode=0



    jobs=[]

    def my_callback(this_errorcode):#propagates the errorcode, if present to higher level variable

        if this_errorcode != 0:
            print("encountered Error",this_errorcode,", ending this run",flush=True)
            nonlocal errorcode
            errorcode = this_errorcode

    #xy-scanning if no deflections are provided
    if len(deflection_arr_info) == 0:
        for i in range(nx):
            for j in range(ny):

                try:#Parallel execution of scan_z
                    _process= workerpool.apply_async(func=scan_z, args=(thresh,z_ax,sstart,i,j,chg_arr_info,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp,shm_height_arr.name,shm_ini_deflection_array.name,shm_unrelaxed_pos_array.name,shm_error_detected.name),callback=my_callback)
                    jobs.append(_process)
                except:
                    errorcode="Error while starting parallel processes"
                    break
                        
            if errorcode != 0:#errorhandling for pool and shared objects
                workerpool.close()
                workerpool.join()
                shm_height_arr.close()
                shm_height_arr.unlink()
                shm_ini_deflection_array.close()
                shm_ini_deflection_array.unlink()
                shm_unrelaxed_pos_array.close()
                shm_unrelaxed_pos_array.unlink()
                shm_error_detected.close()
                shm_error_detected.unlink()
                return height_arr, errorcode
        
                        
    else:#xy-scanning if deflections are provided
        for i in range(nx):
            for j in range(ny):
    
                try:#parallel execution of scan_z
                    _process= workerpool.apply_async(func=scan_z, args=(thresh,z_ax,sstart,i,j,chg_arr_info,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp,shm_height_arr.name,shm_ini_deflection_array.name,shm_unrelaxed_pos_array.name, shm_error_detected.name), callback=my_callback)
                    jobs.append(_process)
                except:
                    errorcode="Multi fehler"
                    break

            if errorcode != 0:#errorhandling for pool and shared objects
                workerpool.close()
                workerpool.join()
                shm_height_arr.close()
                shm_height_arr.unlink()
                shm_ini_deflection_array.close()
                shm_ini_deflection_array.unlink()
                shm_unrelaxed_pos_array.close()
                shm_unrelaxed_pos_array.unlink()
                shm_error_detected.close()
                shm_error_detected.unlink()
    
                return height_arr, ini_deflection_array, unrelaxed_pos_array, errorcode
                    
            
                    
    #wait until all jobs are completed

    workerpool.close()
    #shutdown the workerpool
    workerpool.join()

    if errorcode !=0:#Errorhandling, close and relaese shared memory objects
        shm_height_arr.close()
        shm_height_arr.unlink()
        shm_ini_deflection_array.close()
        shm_ini_deflection_array.unlink()
        shm_unrelaxed_pos_array.close()
        shm_unrelaxed_pos_array.unlink()
        shm_error_detected.close()
        shm_error_detected.unlink()
        if len(deflection_arr_info)==0:
            return height_arr, errorcode
        else:
            return height_arr, ini_deflection_array, unrelaxed_pos_array, errorcode

    #copy from shared memory to local one, close and release shared memory

    height_arr[:] = height_arr_shm[:]
    ini_deflection_array[:] = ini_deflection_array_shm[:]
    unrelaxed_pos_array[:] = unrelaxed_pos_array_shm[:]


    
    shm_height_arr.close()
    shm_height_arr.unlink()
    shm_ini_deflection_array.close()
    shm_ini_deflection_array.unlink()
    shm_unrelaxed_pos_array.close()
    shm_unrelaxed_pos_array.unlink()
    shm_error_detected.close()
    shm_error_detected.unlink()

    
    
    if len(deflection_arr_info) ==0:
        return height_arr, errorcode
    else:
        return height_arr, ini_deflection_array, unrelaxed_pos_array, errorcode

def print_surface_profile(height_arr, x_ax, y_ax, outfile,ini_deflection_arr, unrelaxed_pos_array): # Writes the Grid parameters and the height array into the outfile

    nx = len(height_arr)
    ny = len(height_arr[0])
    
    f_out=open(outfile,"w")
    f_out.write("%d  %d\n" %(nx,ny))
    f_out.write("%f %f\n%f %f\n" %(x_ax[0],x_ax[1],y_ax[0],y_ax[1]))
    for ix in range(nx):
        for iy in range(ny):
            x = x_ax[0] * ix/nx + y_ax[0] * iy/ny
            y = x_ax[1] * ix/nx + y_ax[1] * iy/ny
            f_out.write("%f  %f  %f\n" %(x, y, height_arr[ix][iy]) )

    f_out.close()

    if len(ini_deflection_arr) !=0:
        f_out=open(outfile + "_final_deflections.txt","w")
        for i, it in enumerate(ini_deflection_arr):
            f_out.write("%f %f %f %f %f %f\n" %(unrelaxed_pos_array[i][0], unrelaxed_pos_array[i][1], unrelaxed_pos_array[i][2],it[0],it[1],it[2]))
        f_out.close()

def read_data(infile): # does a sanity check and reads the data
    f_in = open(infile,"r")
    A = f_in.readlines()
    f_in.close()

    x_vec = [0,0]
    y_vec = [0,0]
    z_ax = 0

    nx=0
    ny=0
    nz=0
    d1_arr=[]

    i=0

    # Moves to the beginning of the data
    try:
        while A[i].find("DATAGRID_3D_") == -1:
            i += 1
    except IndexError:
        print("Begin of the datagrid not found!")
        sys.exit()

    i+= 1
    #Get the Grid discretization
    tp = A[i].split()
    nx = int(tp[0])
    ny = int(tp[1])
    nz = int(tp[2])
    #Get the cell vectors
    i+=1
    offset = [float(x) for x in A[i].split()]
    i+=1
    x_vec = [float(x) for x in A[i].split()]
    i+=1
    y_vec = [float(y) for y in A[i].split()]
    i+=1
    tp = [float(z) for z in A[i].split()]
    #Check if the z-vector is defined correctly
    if tp[0] != 0 or tp[1] != 0:
        print("No surface cell! The surface has to be perpendicular to the z direction, the cell has to be given by the vectors (ax,ay,0),(bx,by,0),(0,0,c)!")
        sys.exit()
    z_ax = tp[2]
    i +=1
    if A[i].strip() == "":
        i+=1
    #Read in the data
    while i < len(A):
        tp = A[i].split()
        for t in tp:
            d1_arr.append(float(t))

        i += 1
        if A[i].find("END_DATAGRID") != -1:
            break


    return nx,ny,nz,x_vec,y_vec,z_ax,d1_arr,offset

def filecheck(infile,pp_dir, outfile):#Checks if all files are present
    if not os.path.isfile(infile):
        print("Charge density infile not found! Exiting...")
        sys.exit()

    outdir="./"
    outname=outfile
    tp = outfile.split("/")
    if len(tp) >1:
        outdir = ""
        for i in range(len(tp)-1):
            outdir += tp[i]
            outdir += "/"
        if not os.path.isdir(outdir):
            print("Outfile directory not found. Exiting...")
            sys.exit()
        outname = tp[-1]
            
    x_file=None
    y_file=None
    z_file=None
    if pp_dir != '':
        x_file=pp_dir + "/PPpos_x.xsf"
        y_file=pp_dir + "/PPpos_y.xsf"
        z_file=pp_dir + "/PPpos_z.xsf"
        if not os.path.isfile(x_file):
            print("Probe particle x deflection file not found! Filename has to be PPpos_x.xsf Exiting...")
            sys.exit()
        if not os.path.isfile(y_file):
            print("Probe particle y deflection file not found! Filename has to be PPpos_y.xsf Exiting...")
            sys.exit()
        if not os.path.isfile(z_file):
            print("Probe particle z deflection file not found! Filename has to be PPpos_z.xsf Exiting...")
            sys.exit()

    return x_file, y_file, z_file, outdir, outname

def create_list_from_input(list_str):#auxillary function for the argparser
    output_list=[]
    if list_str.find(",") != -1:
        list_str.lstrip('[')
        list_str.rstrip(']')
        list_str.strip()
        tp = list_str.split(",")
        output_list=[float(x) for x in tp]

    else:
        list_str.strip()
        tp = list_str.split()
        if len(tp) != 3:
            print("could not recognize range format! provide either a list, or a range")
            sys.exit()
        else:
            start = float(tp[0])
            stop = float(tp[1])
            step = float(tp[2])
            val=start
            if (stop < start and step >0) or (stop > start and step < 0):
                print("List range invalid, stepsize positive but stop negative or vice versa")
                sys.exit()
            if step > 0:
                while val < stop:
                    output_list.append(val)
                    val += step
            else:
                while val > stop:
                    output_list.append(val)
                    val += step
    return output_list

def run(args): #Starts the script 

    #define input parameters
    infile=args.infile_chg
    outfile=args.outfile
    threshold=args.iso_value
    sstart=args.scan_start
    pp_dir=args.pp_deflection_dir
    deflection_offset=args.offset
    R=args.tip_distance
    b_decay=args.b_decay
    tip_flip_cutoff=args.no_tip_flip_cutoff
    Amp = args.Amp / 2.0

    #Sanitycheck
    infile_pp_pos_x, infile_pp_pos_y, infile_pp_pos_z, outdir, outname = filecheck(infile,pp_dir, outfile)

    #more input parameters
    iso_list=[]
    if args.iso_range != '':
        iso_list=create_list_from_input(args.iso_range)

    offset_list=[]
    if args.offset_range != '':
        offset_list=create_list_from_input(args.offset_range)

    decay_list=[]
    if args.decay_range!= '':
        decay_list=create_list_from_input(args.decay_range)
        

    
    print("reading input data",flush=True)
    nx,ny,nz,x_vec,y_vec,z_ax,d1_arr,off_dummy = read_data(infile)
    print("converting density",flush=True)
    d3_arr = d1_to_d3_arr(d1_arr,nx,ny,nz)
    del d1_arr

    print("creating shared memory chgarr")
    shm_chg_arr= shared_memory.SharedMemory(create=True,size=d3_arr.nbytes)
    d3_arr_shm=np.ndarray(d3_arr.shape,dtype=np.single,buffer=shm_chg_arr.buf)
    d3_arr_shm[:] = d3_arr[:]

    del d3_arr
    
    chg_arr_info= [shm_chg_arr.name,nx,ny,nz,x_vec,y_vec,z_ax]
    deflection_arr_info=[]
    deflection_offset_idx=0
    
    if infile_pp_pos_x != None and infile_pp_pos_y != None and infile_pp_pos_z != None:
        print("reading probe-particle tip deflection",flush=True)


        print("reading x component",flush=True)
        pp_nx, pp_ny, pp_nz, pp_x_vec, pp_y_vec,pp_z_vec, pp_d1_arr_x, pp_offset = read_data(infile_pp_pos_x)


        print("reading y component",flush=True)
        pp_nx, pp_ny, pp_nz, pp_x_vec, pp_y_vec,pp_z_vec, pp_d1_arr_y, pp_offset = read_data(infile_pp_pos_y)

        
        print("reading z component",flush=True)
        pp_nx, pp_ny, pp_nz, pp_x_vec, pp_y_vec,pp_z_vec, pp_d1_arr_z, pp_offset = read_data(infile_pp_pos_z)


        print("creating array",flush=True)
        pp_d1_arr = np.ndarray(shape=(len(pp_d1_arr_x),3),dtype=np.single)
        print("Filling array",flush=True)
        for k in range(len(pp_d1_arr_x)):
            pp_d1_arr[k,0] = pp_d1_arr_x[k]
            pp_d1_arr[k,1] = pp_d1_arr_y[k]
            pp_d1_arr[k,2] = pp_d1_arr_z[k]
        print("converting deflections",flush=True)
        del pp_d1_arr_x
        del pp_d1_arr_y
        del pp_d1_arr_z
        pp_d3_arr = d1_to_d3_arr(pp_d1_arr,pp_nx,pp_ny,pp_nz)
        del pp_d1_arr

        print("creating shared memory deflections")
        shm_pp_d3_arr= shared_memory.SharedMemory(create=True,size=pp_d3_arr.nbytes)
        pp_d3_arr_shm= np.ndarray(pp_d3_arr.shape,dtype=np.single,buffer=shm_pp_d3_arr.buf)
        pp_d3_arr_shm[:] = pp_d3_arr[:]

        del pp_d3_arr
        
        deflection_arr_info = [shm_pp_d3_arr.name,pp_nx,pp_ny,pp_nz,pp_x_vec,pp_y_vec,pp_z_vec,pp_offset]

    
        deflection_offset_idx = int(deflection_offset / pp_z_vec * pp_nz)



    ##############################
    # execution for list of variables (iso, decay, offset)
    ###############################
    
    if len(iso_list)==0 and len(offset_list)==0 and len(decay_list)==0:
    
        print("extracting height array",flush=True)
        ini_deflection_arr = []
        unrelaxed_pos_arr = []
        if len(deflection_arr_info) ==0:
            height_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)
            ini_deflection_arr, unrelaxed_pos_arr = [] , []
        else:
            height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)

        if errorcode !=0:
            print("Error! %s. Ending Program" %errorcode)
            shm_chg_arr.close()
            shm_pp_d3_arr.close()
            shm_chg_arr.unlink()
            shm_pp_d3_arr.unlink()
            sys.exit()
        print("printing output %s" %(outfile),flush=True)
        print_surface_profile(height_arr,x_vec,y_vec,outfile,ini_deflection_arr, unrelaxed_pos_arr)
        print("Endcode",errorcode)

    
        if no_iso:
            print("INFO: There are xy-points, where the iso-value has not been surpassed for any z.\n If you are looking at molecules in a vacuum, this is fine. For a molecule on a surface this should not happen, check your iso-value")
        if only_iso:
            print("WARNING: There are xy-points, where the density was not lower than the iso-value for any z. Therefore no height could be extracted.\n This height profile can not be used for the constant iso AFM simulation!!\n Check your isovalue and/or your cell.")

    else:
        os.chdir(outdir)
        if len(offset_list) ==0 and len(decay_list)==0:
            for iso in iso_list:
                if len(deflection_arr_info) ==0:
                    height_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)
                    ini_deflection_arr, unrelaxed_pos_arr = [] , []
                else:
                    height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)

                if errorcode != 0:
                    print("Error at iso %f. %s. Skipping this configuration..." %(iso,errorcode))
                    continue
                    

                output = outname + "_iso_%.9f.txt" %iso
                print("printing output %s" % (output),flush=True)
                print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)


        elif len(offset_list)!=0 and len(decay_list)==0:
            for offset in offset_list:
                os.system("mkdir offset_%f" %offset)
                os.chdir("offset_%f" %(offset))
                print("offset %f" %offset, flush=True)
                deflection_offset_idx = int(offset / pp_z_vec * pp_nz)
                
                if len(iso_list) ==0:
                    if len(deflection_arr_info) ==0:
                        height_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)
                        ini_deflection_arr, unrelaxed_pos_arr = [] , []
                    else:
                        height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)

                    if errorcode != 0:
                        print("Error at offset %f. %s. Skipping this configuration..." %(offset,errorcode))
                        continue
                    
                    output = outname
                    print("printing output %s" %(output),flush=True)
                    print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)


                else:
                    for iso in iso_list:
                        if len(deflection_arr_info) ==0:
                            height_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)
                            ini_deflection_arr, unrelaxed_pos_arr = [] , []
                        else:
                            height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,b_decay,tip_flip_cutoff,Amp)

                        if errorcode != 0:
                            print("Error at offset %f iso %f. %s. Skipping this configuration..." %(offset,iso,errorcode))
                            continue

                        

                        output = outname + "iso_%.9f.txt" %iso
                        print("printing output %s" %(output),flush=True)
                        print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)


                os.chdir("..")

        elif len(offset_list)==0 and len(decay_list)!=0:
            for decay in decay_list:
                os.system("mkdir decay_%f" %decay)
                os.chdir("decay_%f" %decay)
                print("decay %f" %(decay),flush=True)

                if len(iso_list) ==0:
                    if len(deflection_arr_info) ==0:
                        height_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)
                        ini_deflection_arr, unrelaxed_pos_arr = [] , []
                    else:
                        height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)

                    if errorcode != 0:
                        print("Error at decay %f. %s. Skipping this configuration..." %(decay,errorcode))
                        continue

                    output = outname
                    print("printing output %s" %(output),flush=True)
                    print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)


                else:
                    for iso in iso_list:
                        if len(deflection_arr_info) ==0:
                            height_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)
                            ini_deflection_arr, unrelaxed_pos_arr = [] , []
                        else:
                            height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)

                        if errorcode != 0:
                            print("Error at decay %f iso %f. %s. Skipping this configuration..." %(decay,iso,errorcode))
                            continue

                        output = outname + "iso_%.9f.txt" %iso
                        print("printing output %s" %(output),flush=True)
                        print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)

                os.chdir("..")

        else:
            for offset in offset_list:
                os.system("mkdir offset_%f" %offset)
                os.chdir("offset_%f" %offset)
                print("offset %f" %offset, flush=True)
                deflection_offset_idx = int(offset / pp_z_vec * pp_nz)
                for decay in decay_list:
                    os.system("mkdir decay_%f" %decay)
                    os.chdir("decay_%f" %decay)
                    print("decay %f" %(decay),flush=True)

                    if len(iso_list) ==0:
                        if len(deflection_arr_info) ==0:
                            height_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)
                            ini_deflection_arr, unrelaxed_pos_arr = [] , []
                        else:
                            height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, threshold, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)

                        if errorcode != 0:
                            print("Error at offset %f decay %f. %s. Skipping this configuration..." %(offset,decay,errorcode))
                            continue

                            
                        output = outname
                        print("printing output %s" %output,flush=True)
                        print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)


                    else:
                        for iso in iso_list:
                            if len(deflection_arr_info) ==0:
                                height_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)
                                ini_deflection_arr, unrelaxed_pos_arr = [] , []
                            else:
                                height_arr, ini_deflection_arr, unrelaxed_pos_arr, errorcode = compute_height_arr(chg_arr_info, iso, sstart,deflection_arr_info,deflection_offset_idx,R,decay,tip_flip_cutoff,Amp)

                            if errorcode != 0:
                                print("Error at offset %f decay %f iso %f. %s. Skipping this configuration..." %(offset,decay,iso,errorcode))
                                continue

                            output = outname + "iso_%.9f.txt" %iso
                            print("printing output %s" %(output),flush=True)
                            print_surface_profile(height_arr,x_vec,y_vec,output ,ini_deflection_arr, unrelaxed_pos_arr)

                    os.chdir("..")

                os.chdir("..")
    shm_chg_arr.close()
    if len(deflection_arr_info) != 0:
        shm_pp_d3_arr.close()
    shm_chg_arr.unlink()
    if len(deflection_arr_info) != 0:
        shm_pp_d3_arr.unlink()
    
            
            
if __name__ == "__main__":
    parser= argparse.ArgumentParser(prog="extract_isosurface_height_map_parallel.py",description="Extracts the height profile of an electron density isosurface. Includes the probe-particle deflections if provided")
    parser.add_argument('infile_chg',help='Charge density in the xsf format')
    parser.add_argument('outfile',help='Name of the height profile outfile')
    parser.add_argument('iso_value',type=float, help='Iso value at which the height will be extracted')
    parser.add_argument('--scan_start','-s',type=float, default=0.8, help='relative height (from 0 to 1) where the height scan begins. Sometimes necessary if there are artifacts in the density')
    parser.add_argument('--pp_deflection_dir','-ppd', type=str, default='',help='Directory where the probe particle deflections are stored. Deflections have to have the filenames PPpos_x.xsf, PPpos_y.xsf and PPpos_z.xsf, as in the ppafm default. If no directory is given, deflections are not included (pure isosurface scan)')
    parser.add_argument('--offset','-o', type=float, default=0.0, help='z-offset between the probe particle deflections and the charge density. If the constant height AFM images are adjusted to fit the charge density, this is not necessary. Can be useful to find the correct parameters for the AFM calculation.')
    parser.add_argument('--tip_distance', '-R', type=float, default=3.0, help='Equilibrium distance between probe-particle and tip. To be chosen as in ppafm')
    parser.add_argument('--b_decay','-b', type=float, default=0.0, help='Decay constant for tunneling current changes due to tip-pp distance changes (shorter distance more current, longer distance less current).')
    parser.add_argument('--no_tip_flip_cutoff','-tfc', action='store_false', help='Disables the detection of unphysical probe particle flips (probe particle higher than tip). This should not be necessary, used only for errorshooting')
    parser.add_argument('--Amp','-A', type=float, default=1.0, help='Peak-to-Peak amplitude over which the charge density is averaged. Should be the same as in ppafm')
    parser.add_argument('--iso_range', '-ir', default='', help='range of iso values to be used. Given either with a , separated list (e.g 0.00002, 0.00003, etc.) or as a range with syntax "start stop stepsize". Overwrites iso_value.')
    parser.add_argument('--offset_range', '-or', default='',  help='range of offsets to be used. Given either with a , separated list (e.g 0.4, 0.5, etc.) or as a range with syntax "start stop stepsize". Overwrites offset.')
    parser.add_argument('--decay_range', '-dr', default='', help='range of tip-pp decays to be used. Given either with a , separated list (e.g 0.05, 0.06, etc.) or as a range with syntax "start stop stepsize". Overwrites b_decay.')

    args=parser.parse_args()
    run(args)
    
