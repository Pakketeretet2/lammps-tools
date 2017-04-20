"""!
\file triangulate.py
\module lammpstools.py

Routines for performing a triangulation and analysing it.
"""



def triangulate( b, rc, dims, method = None ):
    """ ! Triangulates given block. """

    lammpstools = cdll.LoadLibrary("/usr/local/lib/liblammpstools.so")

    pname_base = '/tmp/lammpstools_triangulate_pipe_'
    pname = pname_base + str(os.getpid())

    os.mkfifo(pname)
    triangles = []
    
    try:

        def start_triangulate():
            # Make the child neighborize while the main thread reads from the pipe
            lammpstools.triangulate(void_ptr(b.x), c_longlong(b.meta.N), void_ptr(b.ids),
                                    void_ptr(b.types), c_double(rc),
                                    c_longlong(b.meta.domain.periodic),
                                    void_ptr(b.meta.domain.xlo), void_ptr(b.meta.domain.xhi),
                                    c_longlong(dims), c_longlong(method), c_char_p(pname))
            
        p = Process(target = start_triangulate)
        p.start()

        # Read in the file and store in neighs
        current_l = []
        int_size = 8
        int_fmt  = 'q'
        ints_read = 0
        i = 0
        lc = 1
        with open(pname,"rb") as fifo:
            while True:
                line = fifo.readline()

                if line == '':
                    # Out of data apparently.
                    break
                words = line.split()
                t1 = int(words[0])
                t2 = int(words[1])
                t3 = int(words[2])
                
                triangles.append( [ t1, t2, t3 ] )
                lc += 1
        p.join()

    finally:
        if os.path.exists( pname ):
            os.unlink( pname )
            # time.sleep(0.1)

    return triangles


def triangulation_area( b, triangles ):
    """ Determines the area of given triangulation. """
    im = make_id_map( b.ids )
    A  = 0.0
    it = 1
    for t in triangles:
        # Heron's formula:
        x1 = b.x[t[0]]
        x2 = b.x[t[1]]
        x3 = b.x[t[2]]
        aa = np.linalg.norm( x1 - x2 )
        bb = np.linalg.norm( x1 - x3 )
        cc = np.linalg.norm( x2 - x3 )
        s = (aa + bb + cc) / 2.0
        aaa = math.sqrt( s*(s-aa)*(s-bb)*(s-cc) )
        A  += aaa
        it += 1
    return A
