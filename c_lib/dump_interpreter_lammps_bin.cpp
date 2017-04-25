#include "dump_interpreter_lammps_bin.h"

#include <fstream>


/**
   This interpreter leverages the heavy lifting to a modification of the
   binary2txt tool that comes with LAMMPS.
*/
#include <cstdio>
#include <cstring>


#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#if defined(LAMMPS_SMALLBIG)
typedef int tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#elif defined(LAMMPS_SMALLSMALL)
typedef int tagint;
typedef int bigint;
#define BIGINT_FORMAT "%d"
#else /* LAMMPS_BIGBIG */
typedef int64_t tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#endif



dump_interpreter_lammps_bin::dump_interpreter_lammps_bin(
	const std::string &dname, int file_type )
	: dump_interpreter( dname ), in( nullptr ), id_idx(0), type_idx(1),
	  x_idx(2), y_idx(3), z_idx(4)
{
	in = fopen( dname.c_str(), "rb" );
}

		
dump_interpreter_lammps_bin::~dump_interpreter_lammps_bin()
{
	if( in ) fclose( in );
}


int dump_interpreter_lammps_bin::next_block( block_data &block )
{
	if( !in ) return -1;
	
	if( feof(in) ){
		return 1;
	}
	int size_one, nchunk;
	int status = next_block_meta( block, size_one, nchunk );
	if( !status ){
		return next_block_body( block, size_one, nchunk );
	}else if( status > 1 ){
		// EOF.
		if( status == 1 ) std::cerr << "EOF reached.\n";
		return status;
	}else{
		std::cerr << "Failed to get meta!\n";
		return -1;
	}

}


int dump_interpreter_lammps_bin::next_block_meta( block_data &block,
                                                  int &size_one, int &nchunk )
{
	bigint ntimestep, natoms;
	double xlo[3], xhi[3];
	int boundary[3][2];
	char boundstr[9];
	int triclinic;
	int xy,xz,yz;

	while( true ){
		fread( &ntimestep, sizeof(bigint), 1, in );
		if( !in ) return -1;
		
		if( feof(in) ){
			return 1;
		}

		fread(&natoms,sizeof(bigint),1,in);
		fread(&triclinic,sizeof(int),1,in);
		fread(&boundary[0][0],6*sizeof(int),1,in);
		fread(xlo  ,sizeof(double),1,in);
		fread(xhi  ,sizeof(double),1,in);
		fread(xlo+1,sizeof(double),1,in);
		fread(xhi+1,sizeof(double),1,in);
		fread(xlo+2,sizeof(double),1,in);
		fread(xhi+2,sizeof(double),1,in);
		if (triclinic) {
			fread(&xy,sizeof(double),1,in);
			fread(&xz,sizeof(double),1,in);
			fread(&yz,sizeof(double),1,in);
		}
		fread(&size_one,sizeof(int),1,in);
		fread(&nchunk,sizeof(int),1,in);

		last_meta.N = natoms;
		last_meta.tstep = ntimestep;
		last_meta.xlo[0] = xlo[0];
		last_meta.xlo[1] = xlo[1];
		last_meta.xlo[2] = xlo[2];

		last_meta.xhi[0] = xhi[0];
		last_meta.xhi[1] = xhi[1];
		last_meta.xhi[2] = xhi[2];

		int m = 0;
		for (int idim = 0; idim < 3; idim++) {
			for (int iside = 0; iside < 2; iside++) {
				if (boundary[idim][iside] == 0){
					boundstr[m++] = 'p';
				}else if (boundary[idim][iside] == 1){
					boundstr[m++] = 'f';
				}else if (boundary[idim][iside] == 2){
					boundstr[m++] = 's';
				}else if (boundary[idim][iside] == 3){
					boundstr[m++] = 'm';
				}
			}
			boundstr[m++] = ' ';
		}
		boundstr[8] = '\0';

		last_meta.periodic = 0;
		if( (boundary[0][0] == 0) && (boundary[0][1] == 0) ){
			last_meta.periodic += BIT_X;
		}
		if( (boundary[1][0] == 0) && (boundary[1][1] == 0) ){
			last_meta.periodic += BIT_Y;
		}
		if( (boundary[2][0] == 0) && (boundary[2][1] == 0) ){
			last_meta.periodic += BIT_Z;
		}

		last_meta.boxline  = "ITEM: BOX BOUNDS ";
		last_meta.boxline += boundstr;
		

		return 0;
	}
}

int dump_interpreter_lammps_bin::next_block_body( block_data &block,
                                                  int size_one, int nchunk )
{
	block.copy_meta( last_meta );
	block.resize( block.N );
	
	int maxbuf = 0;
	double *buf = nullptr;
	int n = 0;

	for( int i = 0; i < nchunk; i++ ){
		fread(&n,sizeof(int),1,in);

		// extend buffer to fit chunk size
		
		if( n > maxbuf ){
			if (buf) delete [] buf;
			buf = new double[n];
			maxbuf = n;
		}

		// read chunk and write as size_one values per line

		fread(buf,sizeof(double),n,in);
		n /= size_one;
		block.copy_meta( last_meta );
		block.resize(n);
		int m = 0;
		for( int j = 0; j < n; j++ ){
			for( int k = 0; k < size_one; k++ ){
				// fprintf(fptxt,"%g ",buf[m++]);
				if( k == id_idx ){
					block.ids[j] = buf[m++];
				}else if( k == type_idx ){
					block.types[j] = buf[m++];					
				}else if( k == x_idx ){
					block.x[j][0] = buf[m++];
				}else if( k == y_idx ){
					block.x[j][1] = buf[m++];
				}else if( k == z_idx ){
					block.x[j][2] = buf[m++];
				}
			}
		}
	}
	
	if( buf ) delete [] buf;

	return 0;
}

