#ifndef DUMP_INTERPRETER_LAMMPS_BIN_H
#define DUMP_INTERPRETER_LAMMPS_BIN_H


#include "block_data.h"
#include "dump_interpreter.h"

#include <cstdio>

class dump_interpreter_lammps_bin : public dump_interpreter
{
public:
	dump_interpreter_lammps_bin( const std::string &dname, int file_type = 3 );
	virtual ~dump_interpreter_lammps_bin();
	
	virtual int next_block( block_data &block );

	virtual int next_block_meta( block_data &block, int &size_one, int &nchunk );
	virtual int next_block_body( block_data &block, int size_one, int nchunk );

	virtual bool eof() const
	{
		if ( in ) return std::feof(in);
		else      return false;
	}
	virtual bool good() const
	{
		if( in ) return !std::ferror(in);
		else     return false;
	}
	
private:
	enum { BIT_X = 1,
	       BIT_Y = 2,
	       BIT_Z = 3 };
	
	int id_idx;
	int type_idx;
	int x_idx;
	int y_idx;
	int z_idx;

	int n_cols;
	std::vector<std::string> headers;
	
	std::vector<int>         other_cols;
	std::vector<std::string> other_col_headers;

	int atom_style;
	int scaled;
	std::string last_line;
	block_data last_meta;

	std::FILE *in;
};


#endif // DUMP_INTERPRETER_LAMMPS_BIN_H
