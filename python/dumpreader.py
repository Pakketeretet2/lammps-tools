"""! @package dumpreader
This package contains functions for extracting dump files.

\ingroup lammpstools
"""

import numpy as np
import gzip
import os
import sys
import copy

from dumpreader_cpp import *
from block_data import *



class dumpreader:
    def __init__(self,fname,quiet = True, id_tag = 'id', type_tag = 'type',
                 mol_tag = 'mol'):
        """! Initializes dumpreader with given dump file.

        @param fname Name of dump file to read.
        """
        print("The Python-based dump reader is deprecated and will be removed in a future version. ",
              "Use dumpreader_cpp instead!")

        
        # Check if exists:
        self.at_eof       = 0
        self.bad_file     = 0
        self.n_time_steps = 0
        self.quiet        = quiet
        self.last_block   = None
        self.id_tag       = id_tag
        self.type_tag     = type_tag
        self.x_tag        = None
        self.y_tag        = None
        self.z_tag        = None
        self.fname        = fname

        self.x_tags = [ 'x', 'xs', 'xu' ]
        self.y_tags = [ 'y', 'ys', 'yu' ]
        self.z_tags = [ 'z', 'zs', 'zu' ]
        

        self.scaled_x     = False
        self.scaled_y     = False
        self.scaled_z     = False

        self.molecules    = False
        self.mol_tag      = mol_tag

        if not os.path.isfile(fname):
            print("File %s does not exist!" % fname, file = sys.stderr)
            self.bad_file = 1
            return
        self.smooth_open(fname)

    def __iter__(self):
        return self

    def __next__(self):
        if not self.at_eof:
            return self.getblock()
        else:
            raise StopIteration

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

            if self.x_tag == 'xs': self.scaled_x = True
            if self.y_tag == 'ys': self.scaled_y = True
            if self.z_tag == 'zs': self.scaled_z = True
        


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

            elif w in self.x_tags:
                if self.x_tag is None:
                    self.x_tag = w
                    x_idx = i - 2
                else:
                    if self.x_tag != w:
                        print("x tag changed during dump!", file = sys.stderr)
                        return None
                    else:
                        x_idx = i-2

            elif w in self.y_tags:
                if self.y_tag is None:
                    self.y_tag = w
                    y_idx = i - 2
                else:
                    if self.y_tag != w:
                        print("y tag changed during dump!", file = sys.stderr)
                        return None
                    else:
                        y_idx = i-2
                        
            elif w in self.z_tags:
                if self.z_tag is None:
                    self.z_tag = w
                    z_idx = i - 2
                else:
                    if self.z_tag != w:
                        print("z tag changed during dump!", file = sys.stderr)
                        return None
                    else:
                        z_idx = i-2

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
            
        #print("Forwarded %d lines. Next line is now:" % meta.N, file = sys.stderr)
        #print("'%s'" % line, file = sys.stderr)
        
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

