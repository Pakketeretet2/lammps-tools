"""! @package dumpreader
This package contains functions for extracting dump files.

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy


class domain_data:
    """! A struct containing info about the domain of the dump file"""
    def __init__(self,xlo,xhi,periodic,box_line):
        self.xlo = xlo
        self.xhi = xhi
        self.periodic = periodic
        self.box_line = box_line

class block_meta:
    """! Contains the block "meta"-data, like number of atoms, domain, etc. """
    def __init__(self,t,N,domain):
        self.t = t
        self.N = N
        self.domain = domain
        self.atom_style = "atomic"
    
        
class block_data:
    """! Just a 'struct' that contains some basic data."""
    def __init__(self,meta,ids,types,x, mol = None):
        self.meta = meta
        
        self.ids = ids
        self.types = types
        self.x = x
        if mol is None:
            self.mol = None # np.zeros( meta.N, dtype = int )
        else:
            self.mol = mol
            
        self.other_cols = []


class dump_col:
    """! For if you have/want additional columns in your dump file."""
    def __init__(self,header,data):
        self.header = header
        self.data   = data
        self.N      = len(data)

    


def append_col(block, col):
    """! Appends column to other_cols of block"""
    if( block.meta.N != col.N):
        print("Size mismatch between column and block info!")
        return
    block.other_cols.append(col)
    return block

def read_data(fname):
    """! Reads a data file, converts it to block_data and returns that. """
    with open(fname) as fin:
        
        N = 0
        t = 0
        xlo = np.array([0,0,0],dtype=np.float64)
        xhi = np.array([0,0,0],dtype=np.float64)
        dom = domain_data(xlo, xhi, 0, "")

        skip_two = 0
        
        while True:
            line = fin.readline()
            if skip_two < 2:
                skip_two += 1
                continue
            
            if len(line) == 0:
                # At EOF.
                break
        
            line = line.rstrip()
            if not len(line) == 0:
                # Split line in words:
                words = line.split()
                if( len(words) > 1 ):
                    if words[1] == "atoms":
                        N = int(words[0])

                        ids = np.zeros(N,dtype=int)
                        types = np.zeros(N,dtype=int)
                        x     = np.zeros([N,3], dtype=float)
        
                if( len(words) > 3 ):
                    if( words[2] == "xlo" and words[3] == "xhi" ):
                        dom.xlo[0] = float(words[0])
                        dom.xhi[0] = float(words[1])
                    if( words[2] == "ylo" and words[3] == "yhi" ):
                        dom.xlo[1] = float(words[0])
                        dom.xhi[1] = float(words[1])
                    if( words[2] == "zlo" and words[3] == "zhi" ):
                        dom.xlo[2] = float(words[0])
                        dom.xhi[2] = float(words[1])
                if( words[0] == "Atoms" ):
                    # Read the entire Atoms block.
                    line = fin.readline()
                    jjj = 0
                    
                    while jjj < N:
                        line = fin.readline()
                        line = line.rstrip()

                        words = line.split()
                        idi = int(words[0])
                        ity = int(words[1])
                        xi  = float(words[2])
                        yi  = float(words[3])
                        zi  = float(words[4])

                        ids[jjj] = idi;
                        types[jjj] = ity;
                        x[jjj,0] = xi
                        x[jjj,1] = yi
                        x[jjj,2] = zi
                        
                        jjj += 1

    bm = block_meta(t,N,dom)
    b = block_data(bm,ids=ids,types=types,x=x)
    return b



def read_raw_pts(fname):
    """! Reads a raw text file, assuming coords are x y z in first three columns. """
    
    xyz = np.loadtxt(fname)
    N = len(xyz[:,0])
    t = 0

    xmin = np.min(xyz[:,0])
    xmax = np.max(xyz[:,0])
    ymin = np.min(xyz[:,1])
    ymax = np.max(xyz[:,1])
    zmin = np.min(xyz[:,2])
    zmax = np.max(xyz[:,2])

    ids = np.zeros(N)
    for i in range(0,N):
        ids[i] = i+1
    types = np.ones(N)
    
    xlo = np.array([xmin-2,ymin-2,zmin-2],dtype=np.float64)
    xhi = np.array([xmax+2,ymax+2,zmax+2],dtype=np.float64)
    dom = domain_data(xlo, xhi, 0, "")

    b = block_data(t=0,N=N,domain=dom,ids=ids,types=types,x=xyz)
    return b



class dumpreader:
    def __init__(self,fname,quiet = True, id_tag = 'id', type_tag = 'type',
                 x_tag = 'x', y_tag = 'y', z_tag = 'z', mol_tag = 'mol'):
        """! Initializes dumpreader with given dump file.

        @param fname Name of dump file to read.
        """
        # Check if exists:
        self.at_eof       = 0
        self.bad_file     = 0
        self.n_time_steps = 0
        self.quiet        = quiet
        self.last_block   = None
        self.id_tag       = id_tag
        self.type_tag     = type_tag
        self.x_tag        = x_tag
        self.y_tag        = y_tag
        self.z_tag        = z_tag
        self.fname        = fname

        self.scaled_x     = False
        self.scaled_y     = False
        self.scaled_z     = False

        self.molecules    = False
        self.mol_tag      = mol_tag

        if x_tag == 'xs': self.scaled_x = True
        if y_tag == 'ys': self.scaled_y = True
        if z_tag == 'zs': self.scaled_z = True
        
        if not os.path.isfile(fname):
            print("File %s does not exist!" % fname, file = sys.stderr)
            self.bad_file = 1
            return
        self.smooth_open(fname)
        

    def smooth_open(self,fname):
        """! Gracefully opens file, checks if gzipped or not. """
        # Check if zipped:
        if fname[len(fname)-2:] == "gz":
            print("Opening zipped file %s" % fname, file = sys.stderr)
            self.dump_file = gzip.GzipFile(fname)
        else:
            print( "Opening plain file %s" % fname, file = sys.stderr)
            self.dump_file = open(fname)

    def get_block_meta(self):
        line = self.dump_file.readline()
        
        if len(line) == 0:
            at_eof = 1;
            print("Line was empty! Returning None!", file = sys.stderr)
            return None, None;
        line = line.rstrip()

        t = N = -1
        dom = domain_data( np.array( [ 0, 0, 0 ], dtype=float),
                           np.array( [ 0, 0, 0 ], dtype=float), 0, "pp pp pp" )

        while True:
            if line == "ITEM: TIMESTEP":
                line = self.dump_file.readline().rstrip()
                t = int(line)
                if not self.quiet:
                    print("At t = %d" % t, file = sys.stderr, end = "")
            elif line == "ITEM: NUMBER OF ATOMS":
                line = self.dump_file.readline().rstrip()
                N = int(line)
                if not self.quiet:
                    print(", got %d atoms" % N, file = sys.stderr, end = "")
            elif line.startswith( "ITEM: BOX BOUNDS " ):
                # Extract the periodicity info
                dom.box_line = line
                box_info = line.split()
                if   box_info[3] == "pp": dom.periodic += 1
                if   box_info[4] == "pp": dom.periodic += 2
                if   box_info[5] == "pp": dom.periodic += 4
            
                line = self.dump_file.readline().rstrip()
                data = line.split()
                dom.xlo[0] = float(data[0])
                dom.xhi[0] = float(data[1])
                
                line = self.dump_file.readline().rstrip()
                data = line.split()
                dom.xlo[1] = float(data[0])
                dom.xhi[1] = float(data[1])
                
                line = self.dump_file.readline().rstrip()
                data = line.split()
                dom.xlo[2] = float(data[0])
                dom.xhi[2] = float(data[1])
                if not self.quiet:
                    print(", box = [%1.2f,%1.2f] x [%1.2f,%1.2f] x [%1.2f,%1.2f], periodic = %d" % \
                        ( dom.xlo[0], dom.xhi[0], dom.xlo[1], dom.xhi[1],
                          dom.xlo[2], dom.xhi[2], dom.periodic ),
                          file = sys.stderr, end = "")

            elif line.startswith("ITEM: ATOMS"):
                #if not self.quiet:
                meta = block_meta( t, N, dom )
                return meta, line
            else:

                print("line '", line, "' not recognized! Returning None!", file = sys.stderr )
                return None, None

            line = self.dump_file.readline().rstrip()
        return None, None
            
            
    def getblock_impl(self):
        """! \private Returns the next block in the dump file. """
        
        # Check if at EOF:
        xlo = np.array([0,0,0],dtype=np.float64)
        xhi = np.array([0,0,0],dtype=np.float64)
        dom = domain_data(xlo, xhi, 0, "")
        
        if self.at_eof:
            print("At EOF already! Returning None!", file = sys.stderr)
            return None

        meta, line = self.get_block_meta()
        if meta is None or line is None:
            return None

        xlo[0] = meta.domain.xlo[0]
        xhi[0] = meta.domain.xhi[0]
        xlo[1] = meta.domain.xlo[1]
        xhi[1] = meta.domain.xhi[1]
        xlo[2] = meta.domain.xlo[2]
        xhi[2] = meta.domain.xhi[2]
        
        # Find the right columns
        words = line.split(' ')
        i_idx = t_idx = x_idx = y_idx = z_idx = -1
        mol_idx         = -1
        other_col_idx   = []
        other_col_heads = []
        if not self.quiet:
            print("I got ", word_count, " words.")
            print("They are:", end="")
            for w in words:
                print(" ", w, end = "")
            print("")

        for w, i in zip(words, range(0,len(words))):
            
            if w == self.id_tag:
                i_idx = i - 2
            elif w == self.type_tag:
                t_idx = i - 2
            elif w == self.x_tag:
                x_idx = i - 2
            elif w == self.y_tag:
                y_idx = i - 2
            elif w == self.z_tag:
                z_idx = i - 2
            elif w == self.mol_tag:
                mol_idx = i - 2
                self.molecules = True

            elif( i >= 2 ):
                other_col_heads.append(w)
                other_col_idx.append(i-2)


        if any( x for x in [ i_idx, t_idx, x_idx, y_idx, z_idx ] if x < 0 ):
            print(i_idx, " ", t_idx, " ", x_idx, " ", y_idx, " ", z_idx)
            print("Some index not set! Aborting!")
            return None

        x          = np.zeros( [meta.N, 3], dtype = float )
        ids        = np.zeros( meta.N, dtype = int )
        types      = np.zeros( meta.N, dtype = int )
        n_other_cols = len(other_col_heads)

        if self.molecules:
            mol = np.zeros( meta.N, dtype = int )


        #if n_other_cols > 0:
        #    print >> sys.stderr, "Got %d other cols. They are: " % n_other_cols,
        #    for ch in other_col_heads:
        #        print >> sys.stderr, "%s" % ch,
        #    
        #    print >> sys.stderr, ""

        other_cols = []
        for jj in range(0,n_other_cols):
            other_cols.append( np.zeros(meta.N) )

        # Construct a reverse-lookup array, much faster than ifs.
        lookup_core  = np.zeros(5, dtype=int) # Looks up id, type, x, y, z
        lookup_other = np.zeros(n_other_cols, dtype=int)

        lookup_core[0] = i_idx
        lookup_core[1] = t_idx
        lookup_core[2] = x_idx
        lookup_core[3] = y_idx
        lookup_core[4] = z_idx

        lookup_mol     = mol_idx


        for jj, o_idx in zip(range(0,n_other_cols),other_col_idx):
            lookup_other[jj] = o_idx
        
        for i in range(0,meta.N):
            line = self.dump_file.readline().rstrip()
            words = line.split()

            ids[i]   = int( words[ lookup_core[0] ] )
            types[i] = int( words[ lookup_core[1] ] )
            x[i][0]  = float( words[ lookup_core[2] ] )
            x[i][1]  = float( words[ lookup_core[3] ] )
            x[i][2]  = float( words[ lookup_core[4] ] )

            if self.molecules:
                mol[i] = int( words[ lookup_mol ] )
                if not self.quiet:
                    print("mol[", i, "] = ", mol[i])
            if self.scaled_x:
                Lx = xhi[0] - xlo[0]
                x[i][0] = Lx * x[i][0] + xlo[0]

            if self.scaled_y:
                Ly = xhi[1] - xlo[1]
                x[i][1] = Ly * x[i][1] + xlo[1]

            if self.scaled_z:
                Lz = xhi[2] - xlo[2]
                x[i][2] = Lx * x[i][2] + xlo[2]

            for jj in range(0, n_other_cols):
                other_cols[jj][i] = float( words[ lookup_other[jj] ] )

        # All done, I think.
        if self.molecules:
            b = block_data( meta, ids, types, x, mol = mol )
        else:
            b = block_data( meta, ids, types, x )
        for jj in range(0, n_other_cols):
            occ = dump_col( other_col_heads[jj], other_cols[jj] )
            b = append_col(b, occ)
        

        return b
        
        

    def getblock(self,N=1,pass_last_if_at_eof=True):
        """! Returns the Nth next block in the dump file.
        @param N The next block to return (1 by default, meaning the next block)
        """

        if self.bad_file:
            print("Called getblock on bad file!", file = sys.stderr)
            return None
        
        if self.at_eof:
            print("At EOF already, cannot get more blocks!", file = sys.stderr)
            return self.last_block
        
        if N <= 0:
            print("I expect N >= 1, returning next block...", file = sys.stderr)
            N = 1
        for i in range(0,N):
            tmp = self.getblock_impl()
            if tmp == None:
                self.at_eof = 1
                if pass_last_if_at_eof:
                    print("Encountered EOF, returning last block (if any)", file = sys.stderr)
                    return self.last_block
                else:
                    print("Encountered EOF, returning None", file = sys.stderr)
                    return None
            else:
                self.last_block = tmp
                b = tmp
        if not self.quiet: print("Returning block at ", b.meta.t, file = sys.stderr)
        return b

    def fast_forward(self,Nblocks=1):
        """! Fast-forwards Nblocks. """
        if self.bad_file:
            print("Called getblock on bad file!", file = sys.stderr)
            return None

        if self.at_eof:
            print("At EOF already, cannot get more blocks!", file = sys.stderr)
            return self.last_block

        counter = 0
        while (not self.at_eof) and (counter < Nblocks):
            meta, line = self.get_block_meta()
            if (meta == None) or (line == None):
                # Gone too far
                print("Forwarded past EOF!", file = sys.stderr)
                self.at_eof = True
                return None
            
            # Now just skip meta.N lines:
            for i in range(0,meta.N):
                line = self.dump_file.readline().rstrip()

            counter += 1
            
        print("Forwarded %d lines. Next line is now:" % meta.N, file = sys.stderr)
        print("'%s'" % line, file = sys.stderr)
        
    def getblock_at_time(self,tstep):
        while not self.at_eof:
            b = self.getblock()
            if b is None: return None
            
            if b.meta.t == tstep: return b
        return None
    
    def get_all_blocks(self):
        """! Returns all blocks in the data file in a list of blocks.
             This can be very memory-heavy!"""
        all_data = []
        while not self.at_eof:
            b = self.getblock(N=1, pass_last_if_at_eof = False )
            if b != None:
                all_data.append(b)
        return all_data

    def reset_dump(self):
        """! Resets the dumpreader to the first line of file. """
        self.smooth_open(self.fname) # Not sure how portable this is...

    def getlastblock(self):
        """! Returns the last block. """
        b = self.getblock()
        while not self.at_eof:
            b = self.getblock()
        return b



def block_to_dump_write(b, fname, write_mode = "w", file_mode = "plain"):
    if file_mode == "plain":
        dump_file = open(fname, write_mode)
    elif file_mode == "gz":
        dump_file = gzip.GzipFile(fname,write_mode)
    else:
        print("File mode %s does not exist!" % file_mode, file = sys.stderr)
        sys.exit(-1)

    print("ITEM: TIMESTEP", file = dump_file)
    print(b.meta.t, file = dump_file)
    print("ITEM: NUMBER OF ATOMS", file = dump_file)
    print(b.meta.N, file = dump_file)
    print(b.meta.domain.box_line, file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[0], b.meta.domain.xhi[0] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[1], b.meta.domain.xhi[1] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[2], b.meta.domain.xhi[2] ), file = dump_file)

    if b.meta.atom_style == "atomic":
        atom_string = "ITEM: ATOMS id type x y z"
    elif b.meta.atom_style == "molecule":
        atom_string = "ITEM: ATOMS id mol type x y z"
    
    for c in b.other_cols:
        if c.N != b.meta.N:
            print("Error! Column size does not match rest of block data!")
        atom_string += " " + c.header
    print(atom_string, file = dump_file)
    for i in range(0,b.meta.N):
        #atom_line = "%d %d %f %f %f" % (b.ids[i], b.types[i],
        #                                b.x[i][0], b.x[i][1], b.x[i][2])
        if b.meta.atom_style == "atomic":
            print(b.ids[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2], end = "", file = dump_file)
        elif b.meta.atom_style == "molecule":
            print(b.ids[i], b.mol[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2], end = "", file = dump_file)


        for c in b.other_cols:
            print(c.data[i], file = dump_file, end = "")
        print("", file = dump_file)
    

    
def block_to_dump(b,fname, write_mode):
    """ !Writes given block data to a dump file. """
    
    if fname[len(fname)-2:] == "gz":
        file_mode = "gz"
    else:
        file_mode = "plain"

    block_to_dump_write(b, fname, write_mode, file_mode)

    
    
def block_to_data(b,fname, overwrite = False):
    """ !Writes given block data to a dump file. """

    write_mode = "w"    

    
    # Determine some important params for the data file:
    ntypes = np.max( b.types )
    natoms = b.meta.N
    
    f = open(fname, write_mode)

    print("LAMMPS data file from lammpstools.py.", file = f )
    print("", file = f )
    print("%d atoms"      % natoms, file = f )
    print("%d atom types" % ntypes, file = f )
    print("", file = f )
    print("%f %f xlo xhi" % (b.meta.domain.xlo[0], b.meta.domain.xhi[0]), file = f )
    print("%f %f ylo yhi" % (b.meta.domain.xlo[1], b.meta.domain.xhi[1]), file = f )
    print("%f %f zlo zhi" % (b.meta.domain.xlo[2], b.meta.domain.xhi[2]), file = f )
    print("", file = f )
    print("Atoms", file = f )
    print("", file = f )

    if( b.meta.atom_style == "molecular" ):
        for i in range(0,b.meta.N):
            print("%d %d %d %f %f %f" % (b.ids[i], b.mol[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2]), file = f )
    else:
        for i in range(0,b.meta.N):
            print("%d %d %f %f %f" % (b.ids[i], b.types[i], b.x[i][0], b.x[i][1], b.x[i][2]), file = f )
    


def block_to_xyz_write(b, fname, write_mode = "w"):
    """! Dumps block to xyz file. """

    fp = open(fname,write_mode)

    print(b.meta.N, file = fp)
    print("Atoms. Timestep ", b.meta.t, file = fp )
    
    # xyz needs to be sorted along ids...

    permutation = np.argsort( b.ids )

    for i in permutation:
        print(b.types[i], b.x[i][0], b.x[i][1], b.x[i][2], file = fp, end = "")

        for c in b.other_cols:
            print(c.data[i], file = fp, end = "")
        print("", file = fp)


def copy_meta_new_block( b_old, Xnew, id_new = None, type_new = None, mol_new = None ):
    """! Makes new block_data object for new positions with old meta. """

    bm   = copy.deepcopy(b_old.meta)
    bm.N = Xnew.shape[0]
    X    = np.zeros( [ bm.N, 3 ], dtype = float )

    ids   = np.ones( bm.N, dtype = int )
    types = np.ones( bm.N, dtype = int )
    mol   = np.ones( bm.N, dtype = int )

    
    for i in range(0,bm.N):
        X[i,:]   = Xnew[i,:]
        
        if id_new is None:
            ids[i] = i+1
        else:
            ids[i] = id_new[i]

        if type_new is None:
            types[i] = 1
        else:
            types[i] = type_new[i]

        if mol_new is None:
            mol[i]   = ids[i]
        else:
            mol[i] = mol_new[i]
    
    bb = block_data( bm, ids, types, X, mol )
    return bb
