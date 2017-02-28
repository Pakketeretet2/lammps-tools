"""! @package dumpreader
This package contains functions for extracting dump local files.

\ingroup lammpstools
"""

# Dump reader for local dumps:
import numpy as np
import gzip
import os
import dumpreader

class block_data_local:
    """! Just a 'struct' that contains some basic data."""
    def __init__(self,t,N,domain,cols):
        self.cols = cols


class dumpreader_local:
    def __init__(self,fname):
        """! Initializes dumpreader with given dump file.

        @param fname Name of dump file to read.
        """
        # Check if exists:
        self.at_eof       = 0
        self.bad_file     = 0
        self.n_time_steps = 0
        
        if not os.path.isfile(fname):
            print("File %s does not exist!" % fname)
            self.bad_file = 1
            return

        # Check if zipped:
        if fname[len(fname)-2:] == "gz":
            print("Opening zipped file %s" % fname)
            self.dump_file = gzip.GzipFile(fname)
        else:
            print("Opening plain file %s" % fname)
            self.dump_file = open(fname)
            

    def getblock_impl(self):
        """! \private Returns the next block in the dump file. """
        
        # Check if at EOF:
        N = 0
        t = 0
        xlo = np.array([0,0,0],dtype=float)
        xhi = np.array([0,0,0],dtype=float)
        dom = dumpreader.domain_data(xlo, xhi, 0, "")
        
        if self.at_eof:
            return

        while not self.at_eof:
            line = self.dump_file.readline()
            
            if len(line) == 0:
                at_eof = 1;
                break;
            line = line.rstrip()

            
            if line == "ITEM: TIMESTEP":
                line = self.dump_file.readline().rstrip()
                t = long(line)
                printf( "At t = %d",t )
            elif line == "ITEM: NUMBER OF ENTRIES":
                line = self.dump_file.readline().rstrip()
                N = long(line)
                printf( ", got %d entries", N )
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
                printf( ", box = [%1.2f,%1.2f] x [%1.2f,%1.2f] x [%1.2f,%1.2f], periodic = %d",
                        dom.xlo[0], dom.xhi[0], dom.xlo[1], dom.xhi[1],
                        dom.xlo[2], dom.xhi[2], dom.periodic )
                
            elif line.startswith("ITEM: ENTRIES"):
                # Count the columns and column headers:
                words = line.rstrip().split()
                other_cols = []
                printf(", headers: ")
                for col_header in words[2:]:
                    dc = dumpreader.dump_col(col_header,np.zeros(N))
                    other_cols.append(dc)
                    printf("%s ", dc.header)
                printf("\n")
                
                for i in range(0,N):
                    line = self.dump_file.readline().rstrip()
                    data = line.split()
                    j = 0
                    for c in other_cols:
                        c.data[i] = float(data[j])
                        j += 1
                block = block_data_local( t = t, N = N, domain = dom, cols = other_cols )
                self.last_block = block
                self.n_time_steps += 1
                return block

            

    def getblock(self,N=1,pass_last_if_at_eof=True):
        """! Returns the Nth next block in the dump file.
        @param N The next block to return (1 by default, meaning the next block)
        """

        if self.bad_file:
            print("Called getblock on bad file!")
            return None
        
        if self.at_eof:
            print("At EOF already, cannot get more blocks!")
            return self.last_block
        
        if N <= 0:
            print("I expect N >= 1, returning next block...")
            N = 1
        for i in range(0,N):
            tmp = self.getblock_impl()
            if tmp == None:
                self.at_eof = 1
                if pass_last_if_at_eof:
                    print("Encountered EOF, returning last block (if any)")
                    return self.last_block
                else:
                    print("Encountered EOF, returning None")
                    return None
            else:
                b = tmp
        return b



    
    def get_all_blocks(self):
        """! Returns all blocks in the data file in a list of blocks.
             This can be very memory-heavy!"""
        all_data = []
        while not self.at_eof:
            b = self.getblock(N=1, pass_last_if_at_eof = False )
            if b != None:
                all_data.append(b)
        return all_data



def block_to_dump_plain(b, fname):
    dump_file = open(fname,"w")

    print("ITEM: TIMESTEP", file = dump_file)
    print(b.meta.t, file = dump_file)
    print("ITEM: NUMBER OF ATOMS", file = dump_file)
    print(b.meta.N, file = dump_file)
    print(b.meta.domain.box_line, file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[0], b.meta.domain.xhi[0] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[1], b.meta.domain.xhi[1] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[2], b.meta.domain.xhi[2] ), file = dump_file)
    atom_string = "ITEM: ATOMS id type x y z"
    for c in b.other_cols:
        if c.N != b.meta.N:
            print("Error! Column size does not match rest of block data!", file = sys.stderr)
        atom_string += " " + c.header
    print(atom_string, file = dump_file)
    for i in range(0,b.meta.N):
        atom_line = "%d %d %f %f %f" % (b.ids[i], b.types[i],
                                        b.x[i][0], b.x[i][1], b.x[i][2])
        for c in b.other_cols:
            atom_line += " %f" % c.data[i]
        print(atom_line, file = dump_file)
    

def block_to_dump_gzip(b, fname):
    dump_file = gzip.GzipFile(fname,"w")

    print("ITEM: TIMESTEP", file = dump_file)
    print(b.meta.t, file = dump_file)
    print("ITEM: NUMBER OF ATOMS", file = dump_file)
    print(b.meta.N, file = dump_file)
    print(b.meta.domain.box_line, file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[0], b.meta.domain.xhi[0] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[1], b.meta.domain.xhi[1] ), file = dump_file)
    print("%f  %f" % ( b.meta.domain.xlo[2], b.meta.domain.xhi[2] ), file = dump_file)
    atom_string = "ITEM: ATOMS id type x y z"
    for c in b.other_cols:
        if c.N != b.meta.N:
            print("Error! Column size does not match rest of block data!", file = sys.stderr)
        atom_string += " " + c.header
    print(atom_string, file = dump_file)
    for i in range(0,b.meta.N):
        atom_line = "%d %d %f %f %f" % (b.ids[i], b.types[i],
                                        b.x[i][0], b.x[i][1], b.x[i][2])
        for c in b.other_cols:
            atom_line += " %f" % c.data[i]
        print(atom_line, file = dump_file)

    
def block_to_dump(b,fname, overwrite = False):
    """ !Writes given block data to a dump file. """
    
    if os.path.exists(fname):
        if os.path.isfile(fname):
            if not overwrite:
                print("File already exists! Not overwriting!")
                return
        else:
            print("Given file name is already a dir! Not writing!")
            return
    
    if fname[len(fname)-2:] == "gz":
        block_to_dump_gzip(b,fname)
    else:
        block_to_dump_plain(b,fname)

    
    
