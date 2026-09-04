import sys
import argparse

def idx(ind,ny,nx):#gets index for 1d to 3d array conversion
    iz = ind //(ny*nx)
    rst = ind % (ny*nx)
    iy = rst // nx
    ix = rst % nx

    return ix, iy, iz

def read_in_3d(infile_grid):#Reads in the constant height AFM frequency shifts

    f_in=open(infile_grid,"r")
    A=f_in.readlines()
    f_in.close()

    i=0

    while A[i].find("BEGIN_DATAGRID_3D") == -1:
        i += 1

    i+=1
    tp = [int(x) for x in A[i].split()]

    nx = tp[0]
    ny = tp[1]
    nz = tp[2]

    i+=1
    v0 = [float(x) for x in A[i].split()]
    vs=[]
    for tp in range(3):
        i+=1
        vs.append([float(x) for x in A[i].split()])

    grid=[]
    for j in range(nx):
        grid.append([])
        for k in range(ny):
            grid[-1].append([])
            for l in range(nz):
                grid[-1][-1].append(0)

    for k in range(nx*ny*nz):
        i += 1
        ix, iy, iz = idx(k,ny,nx)
        grid[ix][iy][iz] = float(A[i].split()[0])

    print(nx,ny,nz)
    print(grid[0][0][3])
        
    return grid, v0, vs, nx,ny,nz

def idx_2d(i,ny):#gets index for 1d to 2d array conversion
    ix = i // ny
    iy = i % ny

    return ix, iy

def read_in_height(infile,offset):#Reads in height from the cc-STM height profile
    f_in=open(infile,"r")
    A=f_in.readlines()
    f_in.close()

    i=0

    tp=A[i].split()
    nx = int(tp[0])
    ny = int(tp[1])

    vecs=[]
    for k in range(2):
        i +=1
        vecs.append([float(x) for x in A[i].split()])

    h_arr=[]
    for j in range(nx):
        h_arr.append([])
        for k in range(ny):
            h_arr[-1].append(0)

    max_h=0
    min_h=1000
            
    for k in range(nx*ny):
        ix , iy = idx_2d(k,ny)
        i += 1
        h_arr[ix][iy] = float(A[i].split()[2]) + offset
        if h_arr[ix][iy] < min_h:
            min_h=h_arr[ix][iy]
        if h_arr[ix][iy] > max_h:
            max_h = h_arr[ix][iy]

    oberfl_thr = min_h + (max_h-min_h)/4

    return h_arr, nx, ny, vecs, oberfl_thr


def linearzerlegung(pos,vecs,nx,ny,debug=False):#Transforms coordinates to the height profile grid
    mat = [[vecs[0][0],vecs[1][0]],[vecs[0][1],vecs[1][1]]]
    #print(mat)
    #print("Input_pos:",pos)
    wrk_pos = pos
    if vecs[0][0] == 0.0 or vecs[1][1] == 0.0:
        tp_mat = [mat[1],mat[0]]
        mat = tp_mat
        tp_pos = [pos[1],pos[0]]
        wrk_pos = tp_pos

    if debug:
        print("mat", mat)
        print("wrk_pos",wrk_pos)

    fac = -mat[1][0] / mat[0][0]
    #mat[1][0] += mat[0][0]*fac
    mat[1][1] += mat[0][1]*fac
    wrk_pos[1] += wrk_pos[0]*fac
    if debug:
        print("fac",fac)
        print("mat",mat)
        print("wrk_pos",wrk_pos)

    #print ("fac:",fac)
    
    fac = -mat[0][1] / mat[1][1]
    #mat[0][1] += mat[1][1]*fac
    #mat[0][0] += mat[1][0]*fac
    wrk_pos[0] += wrk_pos[1]*fac

    if debug:
        print("fac",fac)
        print("mat",mat)
        print("wrk_pos",wrk_pos)
    
    #print ("fac",fac)

    #print(wrk_pos[1])
    wrk_pos[0] /= mat[0][0]
    wrk_pos[1] /= mat[1][1]

    #print(wrk_pos[1])
    if debug:
        print("wrk_pos",wrk_pos)
    
    idx_x = wrk_pos[0]*nx
    idx_y = wrk_pos[1]*ny

    if debug:
        print("idx_x",idx_x)
        print("idx_y",idx_y)

    return idx_x, idx_y

    

def get_height_profile_index(ix,iy,nx,ny,nz,v0,vec_afm,h_grid,h_nx,h_ny,h_vec,no_height_interpolation): #Finds the height and the fitting z index for the AFM grid for a given xy AFM index
    pos=[ix/nx * vec_afm[0][i] + iy/ny * vec_afm[1][i] + v0[i] for i in range(2)]#get coordinate (AFM grid)
    #print("Hier ist pos",pos)
    #pos_save= pos
    errorcode=0
    h_ix, h_iy = linearzerlegung(pos,h_vec,h_nx,h_ny) #get height profile indices
    d_x = h_ix -int(h_ix) #difference to next points in height profile, prepare for interpolating
    d_y = h_iy - int(h_iy)
    if no_height_interpolation:
        d_x=0
        d_y=0
    h_ix = int(h_ix)
    h_iy = int(h_iy)
    

    if h_ix < 0 or h_iy < 0:
        print("Index error at height indices h_ix, h_iy", h_ix, h_iy)
        print("Corresponding to Pos: ",pos)
        print("Check if the xy-scan range of the AFM is smaller (or equal) to the area of the charge density cell!")
        #sys.exit()
        return 0,0,0, "index error xy"
    
        
    try:#Interpolate height
        if d_x !=0:
            dh_x = (h_grid[h_ix +1][h_iy] - h_grid[h_ix][h_iy]) * d_x
        else:
            dh_x = 0
        if d_y !=0:
            dh_y = (h_grid[h_ix][h_iy + 1] - h_grid[h_ix][h_iy]) * d_y
        else:
            dh_y = 0
        height = h_grid[h_ix][h_iy] - v0[2] + dh_x + dh_y
    except IndexError:
        print("Index error at height indices h_ix, h_iy", h_ix, h_iy)
        print("Corresponding to Pos: ",pos)
        print("height field infos: nx", h_nx, "ny", h_ny, "vecs", h_vec)
        print("Afm field infos: nx", nx, "ny", ny, "nz", nz, "vecs", vec_afm)
        h_ix, h_iy = linearzerlegung(pos,h_vec,h_nx,h_ny,True)
        print("Check if the xy-scan range of the AFM is smaller (or equal) to the area of the charge density cell!")
        #sys.exit()
        return 0,0,0, "index error xy"

    
    afm_height_idx = height / vec_afm[2][2] * nz #gets the index for the AFM frequency field
    dz = afm_height_idx - int(afm_height_idx)  #Gets difference to integer index for interpolating
    afm_height_idx = int(afm_height_idx)
    height = afm_height_idx * vec_afm[2][2] / nz + v0[2]

    return afm_height_idx, height, dz, errorcode #returns integer index, difference and height at integer index


def get_frequency_field(afm_grid,nx,ny,nz,v0,vec_afm,h_grid,h_nx,h_ny,h_vec,oberfl_thr,no_freq_interpolation,no_height_interpolation):#Gets the frequency shifts for the constant current AFM
    #initializes arrays
    freq_field=[]
    height_field=[]
    hidx_field=[]
    for  j in range(nx):
        freq_field.append([])
        height_field.append([])
        hidx_field.append([])
        for k in range(ny):
            freq_field[-1].append(0)
            height_field[-1].append(0)
            hidx_field[-1].append(0)

    plot_min=10000
    plot_max=-1111
            
    for j in range(nx):#Scanning starts
        for k in range(ny):
            idx,height, dz,errorcode = get_height_profile_index(j,k,nx,ny,nz,v0,vec_afm,h_grid,h_nx,h_ny,h_vec,no_height_interpolation) #Gets z index from height
            if errorcode !=0:#Errorhandling
                return [],[],[],0,0,errorcode
            
            if height < -99: # No iso surface present!                                                                                                                                                                                                   
                freq_field[j][k] = 0
            elif height > 99: # charge density always larger than iso value
                print("Error: The charge density is always larger than the iso value for a given z column, the height profile can not be used.")
                #sys.exit()
                return [],[],[],0,0,"min iso error"

            else:
                
                try:#Interpolate frequency
                    df = (afm_grid[j][k][idx + 1] - afm_grid[j][k][idx]) * dz
                    if no_freq_interpolation:
                        df = 0
                    freq_field[j][k] = afm_grid[j][k][idx] + df
                    height_field[j][k] = height
                    hidx_field[j][k] = idx
                except IndexError:
                    print("Index error at j k idx", j,k,idx)
                    print("Check if the height",height,"is covered in the AFM grid! If not, check scan parameters!")
                    #sys.exit()
                    return [],[],[],0,0,"index z error"


            if height > oberfl_thr:
                #print(j,k)
                if freq_field[j][k] < plot_min:
                    plot_min = freq_field[j][k]
                if freq_field[j][k] > plot_max:
                    plot_max = freq_field[j][k]
                
    return freq_field, height_field, hidx_field, plot_min, plot_max, errorcode

def transpose_matrix(mat):#Transpose matrix for plotting purposes
    nx= len(mat)
    ny = len(mat[0])
    mat_new = []
    for i in range(ny):
        mat_new.append([])
        for k in range(nx):
            mat_new[-1].append(mat[k][i])
        

    return mat_new

def plot_field(outfile,freq_field,plot_min,plot_max, extent):#Plot height or frequency field
    import matplotlib.pyplot as plt

    field = transpose_matrix(freq_field)
    plt.figure(figsize = (8,8))
    #plt.imshow(freq_field,origin="lower",cmap="gray",vmin=plot_min,vmax=plot_max)
    #plt.imshow(field,origin="lower",interpolation="bicubic",cmap="gray",vmin=plot_min,vmax=plot_max)
    plt.imshow(field,origin="lower",interpolation="bicubic",cmap="gray", extent=extent)
    plt.colorbar()
    plt.savefig(outfile,bbox_inches="tight")
    plt.close()

def write_afm_as_txt(freq_field,outfile,vecs_afm,nx_afm,ny_afm):#Writes out frequency field as txt
    f_out=open(outfile,"w")
    for i in range(len(freq_field)):
        for j in range(len(freq_field[0])):
            x = vecs_afm[0][0]*i/nx_afm + vecs_afm[1][0] * j/ny_afm
            y = vecs_afm[0][1]*i/nx_afm + vecs_afm[1][1] * j/ny_afm
            f_out.write("%f %f %f\n" %(x,y,freq_field[i][j]))
            
    f_out.close()
    
def plot_effective_height(height_field,name, extent):#Plots the height profile as used in the AFM image
    import matplotlib.pyplot as plt

    field = transpose_matrix(height_field)
    plt.figure(figsize = (8,8))
    plt.imshow(field,origin="lower",cmap="gnuplot", extent=extent)
    plt.colorbar()
    plt.savefig(name,bbox_inches="tight")
    plt.close()

def plot_slice(afm_3d_field,slic,outfile):#Plot slice of constant height image (only for testing)
    nx = len(afm_3d_field)                                    
    ny = len(afm_3d_field[0])                                 
                                                              
    pltfield= []                                              
    for j in range(nx):                                       
        pltfield.append([])                                   
        for k in range(ny):                                   
            pltfield[-1].append(afm_3d_field[j][k][slic])     
                                                                  
    plot_effective_height(pltfield,outfile)                              

def print_y_row(field,row,filename): #writes one y row of the AFM image as txt (only for testing)
    f_out=open(filename,"w")
    for x in range(len(field)):
        f_out.write("%d %f\n" %(x,field[x][row]))
    f_out.close()

def print_x_row(field,row,filename): #writes one x row of the AFM image as txt (only for testing)
    f_out=open(filename,"w")
    for y in range(len(field[row])):
        f_out.write("%d %f\n" %(y,field[row][y]))
    f_out.close()

def find_files(fold): # finds all height profiles in the supplied folder
    import glob
    import os
    searchpattern= fold + "/**/*.txt"
    files= glob.glob(searchpattern,recursive=True)
    #print(files)
    i= len(files) -1

    def is_number(lst):
        for it in lst:
            try:
                float(it)
            except ValueError:
                return False
        return True
        
    while i >=0:
        f_in=open(files[i],"r")
        line1=f_in.readline().split()
        line2=f_in.readline().split()
        line3=f_in.readline().split()
        line4=f_in.readline().split()
        f_in.close()
        C_1= len(line1) ==2 and is_number(line1)
        C_2= len(line2) ==2 and is_number(line2)
        C_3= len(line3)	==2 and	is_number(line3)
        C_4= len(line4)	==3 and	is_number(line4)

        if not C_1 or not C_2 or not C_3 or not C_4:
            #print (C_1,C_2,C_3,C_4)
            files.pop(i)
        i -= 1

    print(files)
    
    outfiles=[]
    for it in files:
        ht=os.path.split(it)
        outfiles.append(ht[0] + "/Afm_" + ht[1] + ".png")

    print(outfiles)

    
    return files, outfiles
                    
def run(args):
    #Set arguments
    infile_afm=args.infile_afm
    infile_height=args.infile_height
    outfile = args.outfile
    offset = args.offset

    no_height_interpolation = args.no_height_interpolation
    no_freq_interpolation = args.no_freq_interpolation

    
    
    print('Read in AFM frequency file',flush=True)
    grid_afm, v0, vecs_afm, nx_afm, ny_afm, nz_afm = read_in_3d(infile_afm)

    extent = (v0[0] , vecs_afm[0][0] + vecs_afm[1][0], v0[1], vecs_afm[0][1] + vecs_afm[1][1]) #Achsen fuer Plot
    
    if args.full_folder: #Calculates AFM images for all height profiles in folder
        infiles_height, outfiles= find_files(infile_height)
        for i in range(len(infiles_height)):
            print('Read in height profile file',flush=True)
            grid_height, nx_height, ny_height, vecs_height, oberfl_thr = read_in_height(infiles_height[i],offset)
            extent_height = (0, vecs_height[0][0] + vecs_height[1][0], 0, vecs_height[0][1] + vecs_height[1][1])
            print('Calculate height dependent frequency shifts',flush=True)
            freq_field, height_field, hidx_field, plot_min, plot_max, errorcode = get_frequency_field(grid_afm,nx_afm,ny_afm,nz_afm,v0,vecs_afm,grid_height,nx_height,ny_height,vecs_height,oberfl_thr,no_freq_interpolation,no_height_interpolation)
            if errorcode != 0:
                print("Errorcode",errorcode,"skipping configuration")
                continue
            print('Save output',flush=True)
            plot_field(outfiles[i], freq_field,plot_min,plot_max, extent)
            if args.print_effective_height:
                plot_effective_height(height_field,outfiles[i][:-4]+"_effective_height.png", extent)
            if args.print_effective_index:
                plot_effective_height(hidx_field,outfiles[i][:-4]+"_effective_idx.png", extent)
            if args.no_stm:
                plot_effective_height(grid_height,outfiles[i][:-4]+"_height.png", extent_height)
                if args.stm_as_txt:
                    write_afm_as_txt(height_field,outfiles[i][:-4]+"_STM_height.txt", vecs_afm, nx_afm, ny_afm)

            if args.afm_as_txt:
                write_afm_as_txt(freq_field,outfiles[i][:-4]+".txt", vecs_afm,nx_afm,ny_afm)
                
            #print_y_row(grid_height,180,outfiles[i][:-4]+"_yrow_180_height.txt")
            #print_x_row(grid_height,180,outfiles[i][:-4]+"_xrow_180_height.txt")
    else:#Calculates AFM image for one provided height profile
        print('Read in height profile file')
        grid_height, nx_height, ny_height, vecs_height, oberfl_thr = read_in_height(infile_height,offset)
        extent_height = (0, vecs_height[0][0] + vecs_height[1][0], 0, vecs_height[0][1] + vecs_height[1][1])
        
        print('Calculate height dependent frequency shifts') 
        freq_field, height_field, hidx_field, plot_min, plot_max, errorcode = get_frequency_field(grid_afm,nx_afm,ny_afm,nz_afm,v0,vecs_afm,grid_height,nx_height,ny_height,vecs_height,oberfl_thr,no_freq_interpolation,no_height_interpolation)
        if errorcode !=0:
            print("Errorcode",errorcode,"ending run")
            sys.exit()
        print('Save output')
        #print(plot_min,plot_max)
        plot_field(outfile, freq_field,plot_min,plot_max,extent)
        #plot_slice(grid_afm,60,"Constant_height_60.png")
        if args.print_effective_height:
            plot_effective_height(height_field,outfile[:-4]+"_effective_height.png", extent)
        if args.print_effective_index:
            plot_effective_height(hidx_field,outfile[:-4]+"_effective_idx.png", extent)
        if args.no_stm:
            plot_effective_height(grid_height,outfile[:-4]+"_height.png", extent_height)
            if args.stm_as_txt:
                write_afm_as_txt(height_field,outfile[:-4]+"_STM_height.txt", vecs_afm, nx_afm, ny_afm)
        if args.afm_as_txt:
            write_afm_as_txt(freq_field,outfile[:-4]+".txt", vecs_afm,nx_afm,ny_afm)
        #print_y_row(grid_height,180,outfile[:-4]+"_yrow_180_height.txt")
        #print_x_row(grid_height,180,outfile[:-4]+"_xrow_180_height.txt")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Afm_at_Isosurface.py',description='Generates a constant iso (constant current) AFM image. Requires the AFM frequency shifts at different heights, as well as a height profile')
    parser.add_argument('infile_afm',help='AFM frequency shift file at different heigts, produced with ppafm-plot_results --save_df')
    parser.add_argument('infile_height',help='Height profile file for a specific iso surface, generated with extract_isosurface_height_map.py')
    parser.add_argument('outfile',help='Name for the constant iso AFM image outfile')
    parser.add_argument('--offset', '-o', type=float,default=0.0,help='constant offset in angstrom at which the constant iso AFM image will be calculated. If the constant height AFM images are adjusted to fit the charge density, this is not necessary. Can be useful to find the correct parameters for the AFM calculation.')
    parser.add_argument('--no_stm', action='store_false',help='supresses the output of the corresponding height profile (constant current STM image)')
    parser.add_argument('--print_effective_height', '-peh', action='store_true', help='saves height map used in the AFM image (can differ from image of the --no_stm flag due to the possibly different cells of AFM and STM calculation and due to interpolation) ')
    parser.add_argument('--print_effective_index','-pei', action='store_true', help='same as --print_effective_height, but outputs the used grid indices instead of the actual heights. Used mainly for errorshooting.')
    parser.add_argument('--no_freq_interpolation','-nfi', action='store_true', help='disable frequency interpolation between heights. Used mainly for errorshooting.')
    parser.add_argument('--no_height_interpolation','-nhi', action='store_true', help='disable height interpolation from height profile to AFM grid. Used mainly for errorshooting.')
    parser.add_argument('--full_folder', '-ff', action='store_true',help='searches for all height files a folder (includes subdirectories). Folder is given by infile_height argument. Overwrites outfile')
    parser.add_argument('--afm_as_txt','-aat', action='store_true',help='prints the AFM image as a txt file with the format row column frequency_shift')
    parser.add_argument('--stm_as_txt','-sat', action='store_true',help='prints the STM image as a txt file with the format row column height. Works only if --no_stm is not set')
    
    args=parser.parse_args()
    
    run(args)
