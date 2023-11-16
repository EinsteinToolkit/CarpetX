# A custom plotting function that uses openpmd to extract data
# and generate 2d color plots
import numpy as np
import matplotlib
import os
import re
import argparse
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='plot-zslice', description='Python-Based Plotting Program for Carpet/CarpetX')
parser.add_argument('--out-dir', type=str, default='.', help='Directory to write files to')
parser.add_argument('--data-dir', type=str, default='.', help='Directory to write files to')
parser.add_argument('--first-iteration', type=int, default=0, help='The first iteration to extract')
parser.add_argument('--last-iteration', type=int, default=-1, help='The first iteration to extract')
parser.add_argument('--every', type=int, default=1, help='How frequently to extract iterations')
parser.add_argument('--data-index', type=int, default=0, help='The plane to extract data from')
parser.add_argument('--vmin', type=float, default=None, help='The plot range minimum')
parser.add_argument('--vmax', type=float, default=None, help='The plot range maximum')

parser.add_argument('grid_function')
args = parser.parse_args()

thorn, gf_name = args.grid_function.split("::")
print("THORN:",thorn)
print("GF:",gf_name)
out_dir = args.out_dir
data_dir = args.data_dir
data_index = args.data_index

is_carpet = False
is_carpetx = False
ext = ""
for f in os.listdir(data_dir):
    if f.endswith(".xy.asc"):
        is_carpet = True
    elif f.endswith(".bp"):
        is_carpetx = True
        ext = ".bp"
    elif f.endswith(".bp3"):
        is_carpetx = True
        ext = ".bp3"
    elif f.endswith(".bp4"):
        is_carpetx = True
        ext = ".bp4"
    elif f.endswith(".bp5"):
        is_carpetx = True
        ext = ".bp5"

assert is_carpet ^ is_carpetx, data_dir

def carpetx() -> None:
    global gf_name, thorn
    gf_name = gf_name.lower()
    thorn = thorn.lower()
    basename = None
    for f in os.listdir(data_dir):
        g = re.match(r'(.*)\.it\d+\.bp[345]?$', f)
        if g:
            basename = g.group(1)
            break
    assert basename is not None, "Could not find any appropriate .bp files in data-dir."

    import openpmd_api as io
    fname = f"{data_dir}/{basename}.it%08T" + ext
    print("reading:",fname)
    series = io.Series(fname, io.Access.read_only)
    print("openPMD version: ", series.openPMD)
    if series.contains_attribute("author"):
        print("Author: ",series.author)
    frame = 1
    for index in series.iterations:
        if index < args.first_iteration:
            continue
        if args.last_iteration >= 0 and index > args.last_iteration:
            break
        if (index - args.first_iteration) % args.every != 0:
            continue
        iters = series.iterations[index]
        plt.clf()
        for gf in iters.meshes:
            # Note that we are relying on the data to be in order from coarsest to finest.
            uu = iters.meshes[gf]
            dname = f"{thorn}_{gf_name}"
            if dname not in uu:
                continue
            print("gf: ",gf,", index:",index," -> frame: ",frame,sep="")
            data_raw = uu[dname]
            data = data_raw.load_chunk()

            # Data is listed as z, y, x
            x_index = 2
            y_index = 1
            z_index = 0
            assert uu.axis_labels[z_index] == 'z', str(uu.axis_labels)
            assert uu.axis_labels[y_index] == 'y', str(uu.axis_labels)
            assert uu.axis_labels[x_index] == 'x', str(uu.axis_labels)

            # The z-direction is supposed to be unimportant
            data_index = (data.shape[z_index]-1)//2
            series.flush()

            for chunk in data_raw.available_chunks():

                # Get the grid spacing
                dx = uu.grid_spacing[x_index]
                dy = uu.grid_spacing[y_index]

                # Get the grid offsets
                offset = chunk.offset
                x0 = dx*offset[x_index]
                y0 = dy*offset[y_index]

                # Get the size of the grid
                extent = chunk.extent
                nx = extent[x_index]
                ny = extent[y_index]

                # Create a linear array for each axis
                xv = np.linspace(x0,(nx-1)*dx+x0,nx)
                yv = np.linspace(y0,(ny-1)*dy+y0,ny)

                # Sanity
                assert xv.shape[0] == nx
                assert yv.shape[0] == ny

                # Grab the subset of the data to plot
                pdata = data[data_index][offset[y_index]:offset[y_index]+extent[y_index],
                                         offset[x_index]:offset[x_index]+extent[x_index]]

                # Create an x and y variable suitable for use by pcolor
                x = np.zeros(pdata.shape)
                y = np.zeros(pdata.shape)
                assert pdata.shape[1] == xv.shape[0], f"{pdata.shape[0]} == {xv.shape[0]}"
                assert pdata.shape[0] == yv.shape[0], f"{pdata.shape[1]} == {yv.shape[0]}"
                for i in range(xv.shape[0]):
                    for j in range(yv.shape[0]):
                        x[j,i] = xv[i]
                        y[j,i] = yv[j]

                # Get the minimum and maximum
                vmin = args.vmin
                vmax = args.vmax
                if vmin is None:
                    vmin = np.min(data[data_index])
                if vmax is None:
                    vmax = np.min(data[data_index])

                plt.pcolor(x,y,pdata,vmin=vmin,vmax=vmax)

        # All the overlapping plots should appear in this file.
        fig = os.path.join(out_dir,f"{gf_name}%05d.png" % frame)
        print("Saving fig:",fig)
        plt.savefig(fig)
        frame += 1

def carpet() -> None:
    # # 2D ASCII output created by CarpetIOASCII
    # #
    # 0	0 0 0	0	-1.5 -1.5 0	0.893976418630743
    import matplotlib.pyplot as plt
    import numpy as np
    data = np.genfromtxt(os.path.join(data_dir,f"{gf_name}.xy.asc"))
    print(data.shape)
    time_steps = np.sort(np.unique(data[:,0]))

    x_coord_vals = np.sort(np.unique(data[:,1]))
    y_coord_vals = np.sort(np.unique(data[:,2]))

    x_vals = np.sort(np.unique(data[:,5]))
    y_vals = np.sort(np.unique(data[:,6]))
    print("x_vals:",x_vals)
    print("y_vals:",y_vals)

    vals = data[:,8]

    nx = x_coord_vals.shape[0]
    ny = y_coord_vals.shape[0]

    time_step = time_steps[0]
    out = np.zeros((nx,ny))
    x = np.zeros((nx,ny))
    y = np.zeros((nx,ny))

    time_step = 0

    frame = 0
    vmin = args.vmin
    vmax = args.vmax
    if vmin is None:
        vmin = np.min(data[data_index])
    if vmax is None:
        vmax = np.min(data[data_index])

    def mkframe()->None:
        nonlocal frame
        print("ts:",time_step,"x:",np.min(x),np.max(x),"y:",np.min(y),np.max(y),"out:",np.min(out),np.max(out))
        frame += 1
        plt.pcolor(x,y,out,shading='auto',vmin=vmin,vmax=vmax)
        plt.savefig(os.path.join(out_dir,f"{gf_name}%05d.png" % frame).lower())

    for index in range(data.shape[0]):
        row = data[index,:]
        t = int(row[0])
        if t != time_step:
            if time_step < args.first_iteration:
                pass
            elif args.last_iteration >= 0  and time_step > args.last_iteration:
                pass
            elif (time_step - args.first_iteration) % args.every != 0:
                pass
            else:
                mkframe()
            time_step = t
        i = int(row[1])
        j = int(row[2])
        xv = row[5]
        yv = row[6]
        val = row[8]
        x[i,j] = xv
        y[i,j] = yv
        out[i,j] = val #np.nan_to_num(val, nan=0,posinf=0,neginf=0)
    mkframe()

if is_carpet:
    carpet()
else:
    carpetx()